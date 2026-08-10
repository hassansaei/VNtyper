# Issue 233 Mixed Read Routing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restore lossless BAM/CRAM handling for equal paired reads plus singleton/unpaired reads and add merge-blocking, value-sensitive local and Docker guards that make another success-to-failure oracle inversion visible.

**Architecture:** One atomic production wave changes the complete routing-to-Kestrel vertical interface, so no commit leaves callers broken. In parallel, one lane builds a fail-closed compatibility checker and another makes local/Docker execution use the same typed request, captured result, and validators. After those three green lanes integrate, real BAM/CRAM/FASTQ/adVNTR declarations activate the checker in Make and GitHub Actions, followed by the Docker matrix and final documentation.

**Tech Stack:** Python 3.10+, pytest, samtools, Kestrel 1.0.1/Java 11, JSON, Docker/testcontainers, GitHub Actions, Ruff, mypy, branch coverage, diff-cover, Make.

## Global constraints

- Authorities: GitHub #233, `docs/superpowers/specs/2026-08-10-issue-233-mixed-read-routing-design.md`, current `AGENTS.md`, executable tests, and observed real behavior.
- Work in isolated worktrees based on current `origin/main`; recommended integration branch: `fix/issue-233-mixed-read-routing`.
- Use `superpowers:test-driven-development` for every behavior change. RED means a behavioral assertion fails for the intended reason, never collection, imports, missing fixtures, or environment.
- Use `superpowers:systematic-debugging` for every unexpected failure.
- Every new unit-test file declares `pytestmark = pytest.mark.unit`; unit tests use only `tmp_path` and mocks.
- Preserve Python 3.10, Kestrel 1.0.1, calibrated settings, exact summary step names, cwd propagation, and independent shell quoting.
- Never drop, concatenate, recompress, synthesize, or duplicate a non-empty conversion FASTQ.
- Never weaken parity, genotype assertions, tolerance bounds, markers, timeouts, coverage floors, patch coverage, or Docker tiers.
- Expected values come from fresh final-tree real runs. Do not copy v2.0.3 blindly or widen tolerances to obtain green.
- `pipeline.py` (661 LOC) and `kestrel_genotyping.py` (865 LOC) are over the guideline: extract the touched regions rather than growing them.
- New root `scripts/` code is formatted, linted, mypy-checked, and included in branch coverage.
- Do not create, move, or push a release tag.
- Record each task's base SHA, RED/GREEN commands, commit, review package, and verdict in
  `/tmp/vntyper-issue233-sdd/progress.md`. The repository does not currently ignore a
  `.superpowers/` scratch tree, so never create one and assume it is untracked.

## Frozen interfaces

```python
ReadLayout = Literal["paired", "single", "mixed", "invalid", "empty"]

def route_converted_fastqs(
    produced: tuple[str, str, str, str],
    config: Mapping[str, Any],
) -> tuple[str, ...]: ...  # 1-4 paths, R1/R2/other/single order

def construct_kestrel_command(
    *,
    fastq_files: Sequence[str | Path],
    # existing calibrated arguments unchanged
) -> str: ...

@dataclass(frozen=True)
class KestrelInvocation:
    kmer_size: int
    command: str
    log_file: Path

@dataclass(frozen=True)
class KestrelCommandArguments:
    kestrel_path: str | Path
    reference_vntr: str | Path
    vcf_out: str | Path
    java_path: str
    java_memory: str
    max_align_states: int
    max_hap_states: int
    log_level: str
    sample_name: str
    additional_settings: str

def plan_kestrel_invocations(
    *,
    fastq_files: Sequence[str | Path],
    kmer_sizes: Sequence[int],
    output_dir: Path,
    command_arguments: KestrelCommandArguments,
) -> tuple[KestrelInvocation, ...]: ...

def run_kestrel_stage(
    *,
    fastq_files: tuple[str, ...],
    dirs: Mapping[str, str],
    config: Mapping[str, Any],
    sample_name: str,
    log_level: int,
    cwd: str,
    summary: list[dict[str, Any]],
    summary_file_path: str,
    runner: Callable[..., None],
) -> None: ...
```

`pipeline.py` passes its patchable `run_kestrel` symbol as `runner`; the extracted module does not capture a second unpatchable owner.

Structured routing evidence has one stable log record:

```text
READ_SET_ROUTING {"counts":{"other":0,"r1":14690,"r2":14690,"single":1},"layout":"mixed","selected":["output_R1.fastq.gz","output_R2.fastq.gz","output_single.fastq.gz"]}
```

The JSON is sort-key/canonical, uses basenames rather than environment-specific roots, and is emitted exactly once per conversion routed to Kestrel.

Shared orchestration uses:

```python
@dataclass(frozen=True)
class PipelineRequest:
    input_kind: Literal["bam", "cram", "fastq"]
    input_paths: tuple[Path, ...]
    reference_assembly: str
    output_dir: Path
    threads: int
    log_level: str
    cli_options: tuple[str, ...]
    reference_fasta: Path | None

@dataclass(frozen=True)
class PipelineRunResult:
    exit_code: int
    stdout: str
    stderr: str
```

Every local/Docker runner is `Callable[[PipelineRequest], PipelineRunResult]`.
`build_pipeline_argv(request, path_mapper)` is the sole argv builder: it validates a
positive non-Boolean thread count, a supported log level, input arity, and duplicate or
conflicting options, then maps only paths for the target transport. Shared validators
fingerprint the canonical pre-mapping request and own threads, log level, CLI options,
archive checks, routing-record checks, artifacts, and values; runners only supply a path
mapper and execute. No runner may append a hard-coded thread or logging option.

## Parallel execution map

```text
Task 0 reviewed-doc baseline
  ├─ Task 1 atomic production vertical (routing + Kestrel + pipeline)
  ├─ Task 2 compatibility checker core + historical bootstrap seed
  └─ Task 3 typed shared harness + corrected Docker path mapping
                  │
                  └─ integrate three independently green commits
                       ├─ Task 4 exact Kestrel/CRAM fixture proof
                       └─ Task 5 restore real oracles + activate Make/CI guard
                                      └─ Task 6 Docker matrix/workflow sentinel
                                           └─ Task 7 docs + full verification/review/PR
```

- Tasks 1–3 may be developed concurrently in isolated worktrees; their file sets are
  disjoint and each must pass full `make test-unit` before review.
- Task 1 is atomic only across the production/Python caller interface. Its behavior makes
  the ten stale real-data declarations intentionally RED until Task 5 changes those
  oracles: the nine restored BAM nodes, a5c1 adVNTR, their golden-cohort candidates, and
  Docker quick. Do not run or claim integration/Docker green for the Task-1-only commit;
  integrate Task 5 immediately after Tasks 2–4, and identify this bounded expected-red
  interval in both task reports. No caller-breaking Python commit is allowed.
- Task 2 builds/test-drives the executable core but does not activate it against current expected-failure declarations; Task 5 activates it atomically with restored declarations.
- Task 3 is behavior-preserving and keeps current expected-failure cases green while
  making runner transport identical. It may capture fresh characterization values in
  `/tmp`, but it does not edit `tests/test_data_config.json` or activate stricter live
  declarations; Task 5 owns that atomic transition.
- Tasks 4 and 5 may research/run independently after Tasks 1–3 integrate, but Task 5 owns `tests/test_data_config.json`, Make, and CI activation. Serialize shared-file edits.
- Task 6 exclusively owns Docker test/workflow contracts after Task 5.
- Review every task's `base..head`; never trust an agent report without rerunning its focused and full unit gates.

---

## Task 0: Commit the reviewed specification and plan

**Files:**

- `docs/superpowers/specs/2026-08-10-issue-233-mixed-read-routing-design.md`
- `docs/superpowers/plans/2026-08-10-issue-233-mixed-read-routing-plan.md`

- [ ] Read current complete `AGENTS.md`, issue #233, both artifacts, and all adversarial reports.
- [ ] Resolve all Critical/High review findings and run a scoped exact `claude-opus-5` maximum-effort re-review.
- [ ] Run from the documentation worktree:

  ```bash
  mamba run -n vntyper make test-unit
  mamba run -n vntyper make docs-build
  git diff --check
  ```

  Expected baseline on 2026-08-10: 5,072 unit tests; strict docs build exits 0.

- [ ] Commit before implementation worktrees are created:

  ```text
  docs(issue-233): specify lossless mixed-read routing

  Refs #233
  ```

- [ ] Record this commit as every implementation lane's base. The clean-tree requirement applies after this commit, not while the drafts are untracked.

---

## Task 1: Implement the routing-to-Kestrel vertical atomically

**Files:**

- Modify: `vntyper/scripts/read_layout.py`
- Modify: `vntyper/scripts/pipeline_read_routing.py`
- Modify: `vntyper/scripts/kestrel_command.py`
- Add: `vntyper/scripts/kestrel_execution.py`
- Modify: `vntyper/scripts/kestrel_genotyping.py`
- Add: `vntyper/scripts/pipeline_kestrel.py`
- Modify: `vntyper/scripts/pipeline.py`
- Modify: `tests/support/pipeline_harness.py`
- Modify: `tests/unit/test_read_layout.py`
- Modify: `tests/unit/test_pipeline_read_routing.py`
- Modify: `tests/unit/test_kestrel_filtering.py`
- Modify: `tests/unit/test_kestrel_command_extraction.py`
- Add: `tests/unit/test_kestrel_execution.py`
- Modify: `tests/unit/test_run_kestrel_skip.py`
- Modify: `tests/unit/test_pipeline_read_layout_wiring.py`
- Modify: `tests/unit/test_wave2_pipeline_integration.py`
- Modify: `tests/unit/test_pipeline_cwd.py`
- Add: `tests/unit/test_pipeline_kestrel.py`
- Modify: `tests/unit/test_pipeline_artifact_paths.py`
- Modify: `tests/unit/test_pipeline_guards.py`
- Modify: `tests/unit/test_pipeline_preflight.py`
- Modify: `tests/unit/test_summary_record_step.py`
- Modify: `scripts/golden_cohort/case_expectations.py`
- Modify: `scripts/golden_cohort/matrix.py`
- Modify: `tests/unit/test_golden_cohort_diagnostics.py`

### 1A. Characterize extraction boundaries

- [ ] Before new modules exist, extend existing pipeline-harness tests to pin Kestrel paths, config, sample, log level, cwd, summary name/artifact, stale-VCF behavior, missing-VCF failure, and runner failure. Run them GREEN.
- [ ] Add failing import/symbol tests only after characterization; do not claim a nonexistent module passes on the base.

### 1B. RED: classify invalid separately and route every lossless file

- [ ] Add an exhaustive table in `test_read_layout.py`:

  ```python
  @pytest.mark.parametrize(
      ("counts", "verdict", "selected"),
      [
          ({"r1": 14_690, "r2": 14_690, "other": 0, "single": 1}, "mixed", ("r1.gz", "r2.gz", "single.gz")),
          ({"r1": 3_474, "r2": 3_474, "other": 0, "single": 93}, "mixed", ("r1.gz", "r2.gz", "single.gz")),
          ({"r1": 9, "r2": 9, "other": 4, "single": 2}, "mixed", ("r1.gz", "r2.gz", "other.gz", "single.gz")),
          ({"r1": 0, "r2": 0, "other": 4, "single": 2}, "mixed", ("other.gz", "single.gz")),
      ],
  )
  def test_lossless_layout_routes_every_nonempty_file_once(counts, verdict, selected):
      assert classify_layout(**counts) == verdict
      consumed, stranded = route_fastqs(verdict, PATHS, counts)
      assert consumed == selected
      assert stranded == ()
  ```

- [ ] Add a non-tautological duplicate-path negative control where `other` and `single`
  resolve to the same path; routing must fail rather than double-count one file.
- [ ] Assert `(900,880,0,0)`, `(10,0,0,0)`, and `(0,10,0,0)` are `invalid`, consume none, and strand every non-empty path. Empty and negative counts remain failures/empty as specified.
- [ ] In counted-routing tests, use valid tiny gzip FASTQ records for `(3,3,0,1)`, `(3,3,2,1)`, and `(0,0,2,1)`; assert tuple/order plus the exact canonical routing record. Define “malformed” here narrowly as unreadable/decompression failure or an incomplete four-line group; syntactic FASTQ header/quality validation is out of scope.
- [ ] Run focused tests. RED must be valid mixed rejection and invalid/mixed conflation.

### 1C. RED: N-file Kestrel command and extracted plan

- [ ] Change command tests to `fastq_files`; pin byte-identical one/two-file commands and exact three/four-file commands.
- [ ] Assert exactly one `-sSAMPLE`, ordered operands, no duplicate, no `None`, and independent quoting for spaces, quotes, semicolons, and `$()`.
- [ ] Reject a bare scalar `str` or `Path` container before `tuple(...)`, empty sequence, `None`/empty/unsupported elements, and a relative path whose string begins `-`. Absolute `/dir/-reads.fastq.gz` and relative `dir/-reads.fastq.gz` remain valid. Do not claim shell quoting alone makes `-reads.fastq.gz` positional.
- [ ] Parse `additional_settings` with `shlex.split` and reject sample/input-grouping
  overrides: `-s`, attached `-sNAME`, `--sample`, `--sample=<name>`,
  `--filespersample`, and `--filespersample=<value>`. Otherwise operator config could
  silently split the ordered tuple into multiple samples. Preserve unrelated existing
  additional settings.
- [ ] Implement the exact frozen `KestrelCommandArguments` fields above. `output_dir`
  remains the planner argument used for `output.sam`, `--temploc`, and the per-k log; it
  is not duplicated in the dataclass. Define the exact `KestrelInvocation`/planner
  signature above and reject unsupported field types before command rendering.
- [ ] Test empty/Boolean/non-positive k-mer values, stable order, one invocation per k-mer, immutable input tuple, and exact log path `output_dir / f"kestrel_kmer_{k}.log"`.
- [ ] Run focused tests. RED must be missing N-file/planner behavior.

### 1D. RED: whole-pipeline tuple forwarding

- [ ] Make the harness router return `("R1", "R2", "single")`; assert the same tuple reaches `runner` as `fastq_files` for BAM and CRAM.
- [ ] Add paired/single FASTQ controls proving preprocessing and post-alignment conversion precede routing, then the exact routed tuple reaches Kestrel.
- [ ] Assert route omission, reordering, duplicate file, empty tuple, Kestrel failure, and summary failure do not produce a false successful summary.
- [ ] Pin each conversion summary's FASTQ artifact to `fastq_files[0]`, guarded by the
  non-empty invariant, so tuple migration does not retain a removed `fastq1` variable or
  change the recorded artifact.
- [ ] Capture archive-protected paths from the original CLI inputs before SHARK,
  alignment, or routing reassigns anything. Assert BAM, CRAM, original FASTQ1/FASTQ2,
  reference, BED, and BWA reference reach `create_safe_archive`; the routed tuple itself
  must never be treated as one path.
- [ ] Retarget `scripts/golden_cohort/case_expectations.py`,
  `scripts/golden_cohort/matrix.py`, and `test_golden_cohort_diagnostics.py`: the 32
  measured mixed IDs no longer declare candidate failure and instead use ordinary
  success/artifact expectations. Preserve causal failure checks for a purpose-built
  genuinely `invalid` parity case.
- [ ] Run focused tests. RED must be two-value destructuring/old Kestrel keywords.

### 1E. Minimal atomic GREEN

- [ ] Add `invalid`; route every non-empty valid `mixed` file in R1/R2/other/single order. `route_converted_fastqs` returns `tuple[str, ...]`, rejects stranded/duplicates/non-positive selections, and emits exactly one structured record.
- [ ] Implement typed N-file command/planner and change `run_kestrel` to use it while preserving stale-VCF removal, cwd, command failure, SAM/BAM conversion, result processing, and no-VCF failure.
- [ ] Extract `run_kestrel_stage` with the frozen signature. `pipeline.py` passes `runner=run_kestrel`, keeping the current harness patch owner effective.
- [ ] Leave the existing two-file SHARK, fastp, BWA, and direct-input variables
  untouched. Introduce `kestrel_fastq_files` only at the post-alignment routing boundary
  and replace both conversion destructures, avoiding an unrelated preprocessing refactor.
- [ ] Confirm `pipeline.py` is below 650 LOC or materially reduced by the complete touched stage; confirm `kestrel_genotyping.py` did not grow and the touched pure planning moved out.
- [ ] Preserve import-time `kestrel_config`, the
  `kestrel_genotyping.construct_kestrel_command` re-export identity, and the source-text
  six-gate tripwire. Tests patch the module global as repository policy requires.
- [ ] Run `rg` for every changed symbol and update every caller/test in this atomic
  commit; the focused list may not hide an old keyword or two-value unpack.

### 1F. Verify and commit one green vertical

```bash
mamba run -n vntyper pytest -q \
  tests/unit/test_read_layout.py \
  tests/unit/test_pipeline_read_routing.py \
  tests/unit/test_kestrel_filtering.py \
  tests/unit/test_kestrel_command_extraction.py \
  tests/unit/test_kestrel_execution.py \
  tests/unit/test_run_kestrel_skip.py \
  tests/unit/test_pipeline_kestrel.py \
  tests/unit/test_pipeline_read_layout_wiring.py \
  tests/unit/test_wave2_pipeline_integration.py \
  tests/unit/test_pipeline_cwd.py
mamba run -n vntyper make format-check
mamba run -n vntyper make type-check-all
mamba run -n vntyper make test-unit
```

- [ ] Mutation proof: omit singleton, reorder selected files, accept unequal mates, add a
  second `-s`, allow `--filespersample`, truncate to two Kestrel inputs, pass the tuple
  into archive protection, or retain a mixed-failure golden declaration; each owning
  test fails. Restore and rerun GREEN.
- [ ] Commit only after full unit GREEN:

  ```text
  fix(pipeline): route complete read sets to Kestrel

  Closes #233
  ```

---

## Task 2: Build the compatibility checker and authoritative bootstrap seed

**Files:**

- Add: `scripts/integration_compatibility.py`
- Add: `scripts/check_integration_compatibility.py`
- Add: `tests/compatibility/real_success_baseline.json`
- Add: `tests/unit/test_integration_compatibility.py`
- Add: `tests/unit/test_integration_compatibility_cli.py`

This task builds a green, inactive engine. Task 5 wires it to live declarations and CI atomically.

### 2A. Freeze the contract schema

- [ ] Define identity as `(suite, test_name)`, allowing the same input in BAM and adVNTR suites.
- [ ] Each success contract contains exact:
  - suite/test name;
  - input paths and archive MD5 or derived-CRAM record digest;
  - reference assembly and optional reference path/digest;
  - positive thread count, log level, normalized ordered CLI options, and module set;
  - required routing counts/selected basenames when conversion applies;
  - required artifact paths/archive state;
  - required Kestrel/adVNTR/cross-match field names, expected values, and numeric tolerance no broader than the reviewed declaration;
  - `compatibility_since` and provenance commit.
- [ ] Reject missing/extra keys, duplicate identities, absolute/escaping paths, wrong types, Boolean-as-integer, duplicate options/artifacts, unbounded tolerances, and unresolved resource digests.

### 2B. Bootstrap from authoritative history

- [ ] Pin provenance commit `52c4146596fef2d1e2402991fbab062ba8021889^` and the exact nine BAM plus a5c1 adVNTR identities that were successful there. Tests independently hard-code those ten identities and obtain their historical declarations with `git show`; omission or mutation of any one fails.
- [ ] Add every other qualifying live success (single-end BAM, FASTQ, pure-paired adVNTR) to the manifest. The eventual Task 5 activation enforces bidirectionality: every manifest row is one live success and every real integration success is in the manifest.
- [ ] Bootstrap routing fields from the pinned historical declarations and retained
  design evidence, not from Task 1's parallel implementation output; the literal ten-ID
  seed remains reviewable and independently green.
- [ ] A missing base manifest is accepted only when the current manifest validates against the pinned bootstrap commit and exact ten-ID seed. It is never an unchecked empty bootstrap.

### 2C. Executable checker core

- [ ] Pure core accepts `base_manifest`, `current_manifest`, live test config, and resource config. In one call it:
  1. validates both manifests;
  2. rejects removal or mutation of any base identity;
  3. validates bootstrap history when no base manifest exists;
  4. validates every manifest identity resolves exactly once to an exit-0 live case with matching invocation, artifacts, values, and tolerances;
  5. rejects any qualifying live real success missing from the manifest.
- [ ] The CLI requires `--base-revision` for CI. Local mode may fall back to `git merge-base origin/main HEAD`, printing that it is local-only evidence. It uses `git show BASE:...`; missing/unreachable/all-zero base, shallow history, subprocess failure, malformed JSON, or inconsistent Git output fails closed.
- [ ] Unit tests use temporary Git repositories and cover PR base, branch behind base, direct-push `before` SHA, absent bootstrap manifest, deletion, rename, suite/input/reference/options/modules/digest/artifact/value/tolerance mutation, newly added success without baseline, valid append, and shared input across suites.
- [ ] CLI tests assert argv and exit/error text; every new function/failure path is covered.

### 2D. Verify inactive engine

```bash
mamba run -n vntyper pytest -q \
  tests/unit/test_integration_compatibility.py \
  tests/unit/test_integration_compatibility_cli.py
mamba run -n vntyper ruff check scripts/integration_compatibility.py scripts/check_integration_compatibility.py \
  tests/unit/test_integration_compatibility.py tests/unit/test_integration_compatibility_cli.py
mamba run -n vntyper mypy scripts/integration_compatibility.py scripts/check_integration_compatibility.py
mamba run -n vntyper make test-scripts-cov
mamba run -n vntyper make test-unit
```

- [ ] Commit the inactive, fully tested engine:

  ```text
  test(integration): define append-only success contracts

  Closes #233
  ```

---

## Task 3: Make local and Docker orchestration genuinely identical

**Files:**

- Modify: `tests/support/orchestration.py`
- Modify: `tests/parametrization.py`
- Modify: `tests/integration/test_pipeline_integration.py`
- Modify: `tests/docker/conftest.py`
- Modify: `tests/docker/test_docker_pipeline.py`
- Add: `tests/unit/test_shared_pipeline_orchestration.py`
- Add: `tests/unit/test_docker_pipeline_runner.py`

This task preserves current declared outcomes; it changes transport and validation ownership only.

### 3A. RED: typed request/result and shared ownership

- [ ] Add the frozen `PipelineRequest`/`PipelineRunResult` types. Require local and Docker runners to accept/return them.
- [ ] Shared BAM orchestration constructs the request including exact `threads`,
  `log_level`, and `cli_options`; move archive validation inside shared success
  validation. Add/test `build_pipeline_argv(request, path_mapper)` and assert local and
  Docker receive byte-identical canonical requests before path mapping. Raw host/container
  paths may differ, but normalized resource identities and every non-path argv token must
  match exactly.
- [ ] Add strict shared FASTQ value validation for Kestrel, coverage, summary, and report,
  and strict adVNTR validation for Kestrel, coverage, pipeline summary, exact adVNTR
  fields, cross-match step/result, report, and archive state. Unit-test these validators
  with complete explicit outcome fixtures and record fresh base characterizations in
  `/tmp`; do not add live manifest values or invoke the strict path from mutable live
  declarations in this task.
- [ ] Shared result validation parses exactly one canonical `READ_SET_ROUTING` record for alignment-derived cases and asserts declared counts/selected basenames. Missing, duplicate, malformed, dropped singleton, reordered, or extra selected paths fail.
- [ ] Expected failures declare a stable causal diagnostic. Shared orchestration checks it
  for both local and Docker results; an exit code of 1 is never accepted by itself.
- [ ] Preserve the current declared-failure diagnostic adapter only as a temporary
  compatibility seam, because Task 3 must remain green before the live outcomes change.
  Tests name this seam and prove it is still reached only for current nonzero
  declarations; Task 5 removes it atomically with the restored outcomes.

### 3B. RED: correct Docker path/output mapping

- [ ] Extract/test one host-to-container mapper: every input must resolve under `tests/data`, preserve subdirectories, and map to `/opt/vntyper/input/<relative>`; escapes/symlinks/unregistered paths fail.
- [ ] Output mapper preserves each isolated per-case directory under the mounted `docker_output*` root. It rejects the mount root for parameterized cases, requires an initially empty case directory, and never writes global `/opt/vntyper/output`.
- [ ] Correct the current FASTQ runner, which ignores `output_dir` and maps by basename. Add CRAM request handling with optional explicit reference mapped read-only from either the test-data mount or image path.
- [ ] Capture Docker exec output in `PipelineRunResult` rather than discarding it. Local runner preserves stdout/stderr; shared validator searches their concatenation.
- [ ] Delete local `--threads 4` and Docker `--threads 2` defaults. Thread count and log
  level come only from `PipelineRequest`; mutation tests that reintroduce either
  transport-specific default must fail canonical-request/argv equality.

### 3C. GREEN and verify behavior preservation

```bash
mamba run -n vntyper pytest -q \
  tests/unit/test_shared_pipeline_orchestration.py \
  tests/unit/test_docker_pipeline_runner.py \
  tests/unit/test_integration_outcome_contract.py \
  tests/unit/test_golden_cohort_matrix.py
mamba run -n vntyper make format-check
mamba run -n vntyper make type-check-all
mamba run -n vntyper make test-unit
```

- [ ] Mutation proof: make Docker omit `--fast-mode`, archive options, a nested input subpath, or the case output suffix; each exact request/mapping test fails. Make validator ignore routing output or change FASTQ Confidence; owning tests fail.
- [ ] Commit independently green:

  ```text
  refactor(tests): unify local and Docker pipeline oracles

  Closes #233
  ```

---

## Task 4: Prove exact Kestrel and CRAM fixture behavior

**Files:**

- Add: `tests/integration/test_kestrel_multifile.py`
- Modify: `scripts/make_cram_fixtures.py`
- Modify: `tests/unit/test_integration_tier_hygiene.py`
- Modify: `tests/unit/test_make_cram_fixtures.py`

### 4A. Deterministic fixture generation

- [ ] Run `make cram-fixtures` from registered BAM data. Use:
  - `tests/data/cram/example_b178_hg19_subset.cram`, generated losslessly with `no_ref=1` from `tests/data/example_b178_hg19_subset.bam`;
  - `tests/data/cram/remapped/bwa/hg19/example_b178_hg19_bwa.cram` from the remapped BAM;
  - `tests/data/cram/manifest.json` record count/digest as provenance.
- [ ] Extend the generator with a purpose-built real-read reference-compressed CRAM at
  `tests/data/cram/reference-compressed/example_b178_hg19_subset.cram`. Generate it from
  registered b178 with `samtools view -C -T reference/alignment/chr1.hg19.fa`, index it,
  and prove the explicit-reference decoded record digest equals the source. Manifest
  fields pin encoding=`reference-compressed`, source digest, and reference SHA-256.
- [ ] Unit-test missing/wrong reference, command argv, index/decode failure, digest
  mismatch, manifest fields, and rerun idempotency. Existing `no_ref=1` fixtures remain
  unchanged; they prove CRAM container/routing behavior, not reference decode.
- [ ] In `tests/integration/test_cram_pipeline.py`, execute the generated
  reference-compressed CRAM with the exact reference and assert its decoded-record digest
  equals the source BAM. Execute read-only negative controls with no reference and with a
  wrong tiny reference; both must fail at CRAM decode with a causal diagnostic before the
  pipeline can report success. Mocked generator tests do not satisfy this proof.

### 4B. Exact multi-file Kestrel test

- [ ] In `tmp_path`, call production conversion on registered b178 and 40cf BAMs to produce four FASTQs; do not commit generated data.
- [ ] Invoke the production Kestrel command builder with identical calibrated settings for pair-only and all-nonempty tuples.
- [ ] Parse VCF records as normalized tuples excluding volatile headers. For b178, assert pair-only multiset equals all-file multiset (6,280 observed during design). For 40cf, assert the all-file multiset contains every pair-only record, record the fresh additional count (9 observed), and both processed results remain `Negative`.
- [ ] Assert each command has one sample, exit 0, non-empty VCF/SAM, exact JAR checksum/version, and selected input paths once. No text-only help assertion substitutes for execution.
- [ ] Mark only `integration`; unit hygiene pins marker and registered-data usage.

### 4C. RED/GREEN and commit

- [ ] Before Task 1 integration, the public builder cannot express three files; record that behavioral RED. After Task 1, run:

  ```bash
  VNTYPER_TEST_DATA_SKIP_DOWNLOAD=1 mamba run -n vntyper pytest -q \
    tests/integration/test_kestrel_multifile.py -m integration
  ```

- [ ] Commit:

  ```text
  test(kestrel): prove real multi-file sample inputs

  Closes #233
  ```

---

## Task 5: Restore real success oracles and activate the gate

**Files:**

- Modify: `tests/test_data_config.json`
- Modify: `tests/compatibility/real_success_baseline.json`
- Modify: `tests/parametrization.py`
- Modify: `tests/support/orchestration.py`
- Modify: `tests/integration/test_pipeline_integration.py`
- Add: `tests/integration/test_cram_pipeline.py`
- Modify: `tests/unit/test_integration_outcome_contract.py`
- Modify: `tests/unit/test_integration_compatibility.py`
- Modify: `tests/unit/test_makefile_recipes.py`
- Modify: `tests/unit/test_workflow_consistency.py`
- Modify: `Makefile`
- Modify: `.github/workflows/ci-tests.yml`

### 5A. RED: current bug must fail successful real declarations

- [ ] Restore all nine pre-`52c41465` BAM successes and a5c1 adVNTR success as initial hypotheses, including their previous value/artifact assertions. Add explicit `expected_fastq_records` and `expected_selected_fastqs` to alignment cases.
- [ ] Add current FASTQ characterization values for Kestrel, coverage, summary, and report. Keep b178 paired and 6449 single controls.
- [ ] Atomically switch every successful BAM/CRAM/FASTQ/adVNTR declaration to the strict
  shared outcome schema and validator built in Task 3. Missing value/artifact fields now
  fail rather than taking a compatibility branch. Remove the temporary valid-mixed
  failure adapter; retain a causal diagnostic renderer only for genuinely invalid
  mate-parity cases.
- [ ] Add `integration_tests.cram_tests` with exact keys: `test_name`, `cram`,
  `source_bam`, `fixture_record_digest`, `fixture_reference_sha256` (null only for
  `no_ref=1` fixtures), `reference_assembly`, `reference_fasta`, `threads`, `log_level`,
  `cli_options`, routing counts/selected, artifacts, Kestrel/coverage values, and archive
  state. Add getters/IDs
  and `run_cram_test_case` using `PipelineRequest`; the live/baseline validator binds the
  reference path and digest together.
- [ ] Cases:
  - `b178_hg19_original_derived_mixed`: generated `tests/data/cram/example_b178_hg19_subset.cram`, `(14690,14690,0,1)`;
  - `b178_hg19_reference_compressed_mixed`: generated
    `tests/data/cram/reference-compressed/example_b178_hg19_subset.cram`, exact image
    reference `reference/alignment/chr1.hg19.fa`, decoded-record digest equal to the
    source BAM, `(14690,14690,0,1)`;
  - `b178_hg19_remapped_pure_paired`: `tests/data/cram/remapped/bwa/hg19/example_b178_hg19_bwa.cram`, `(14689,14689,0,0)`.
- [ ] Every real declaration uses `threads: 2` and `log_level: DEBUG`, chosen to match
  the constrained GitHub/Docker runner. Local and Docker runners consume these exact
  fields; neither supplies a different default.
- [ ] Add `example_b178_hg19_subset_default` and
  `example_40cf_hg38_subset_default` without `--fast-mode`, with exact routing, values,
  and artifacts. Existing fast-mode cases remain controls. The default-mode b178 success
  becomes the PR Docker quick sentinel so the reporter's path is protected directly.
- [ ] Run b178 BAM and original-derived CRAM on the base. RED must be exit 1 versus required 0 with exact current refusal, not missing data.

### 5B. GREEN: fresh final-tree local matrix

- [ ] Run all nine BAM cases. Capture exact routing JSON, Kestrel row, coverage, summary, report, archive, and runtime. b178's singleton and 40cf's 93 singletons must appear exactly once in selected inputs.
- [ ] Run all three CRAM cases, paired b178 FASTQ, single 6449 FASTQ, mixed a5c1 adVNTR,
  and pure-paired b178 adVNTR control. Only the reference-compressed case may satisfy
  the image-reference packaging/decode acceptance check.
- [ ] a5c1 asserts Kestrel values, adVNTR columns, cross-match step/result, coverage, report, and archive. Any changed value needs a read-level explanation; never widen tolerance.
- [ ] Update declarations only from captured results and rerun the exact nodes GREEN.
- [ ] Mutation proof: remove one required FASTQ value, adVNTR field, cross-match result,
  archive assertion, or routing count; each exact live case must fail before the
  compatibility baseline is accepted.

### 5C. Activate bidirectional compatibility checking

- [ ] Append every new successful BAM/CRAM/FASTQ/adVNTR identity and full invocation/oracle fingerprint. Checker must now prove baseline→live and qualifying live success→baseline.
- [ ] Add `make check-integration-compatibility`; invoke it from `check-all` and test recipe failure is not masked. The executable performs both live validation and historical append-only comparison in one process.
- [ ] In `.github/workflows/ci-tests.yml`, add a named step in the full-history unit
  matrix with `if: matrix.python-version == '3.12'`, so it runs exactly once per workflow:
  - PR: `--base-revision "origin/${{ github.base_ref }}"`;
  - push: `--base-revision "${{ github.event.before }}"`;
  - fail closed on unsupported/empty/all-zero/unreachable base.
- [ ] Pass the expression through an environment variable and double-quote it in the
  shell; do not interpolate an event value into executable shell text. The step invokes
  `make check-integration-compatibility INTEGRATION_COMPAT_BASE="$BASE_REVISION"`.
- [ ] Add semantic workflow tests pinning event-to-base selection, exact Make invocation, full-history checkout, and absence of `continue-on-error`/`|| true`. Direct-main push and PR tests use real temporary Git history.
- [ ] Because a workflow changes, run `make ci-local` after focused GREEN.

### 5D. Mutation-sensitive verification

- [ ] Independently mutate each of: b178 exit, ID, suite, input digest, reference, fast mode, archive option, module set, routing singleton selection, required artifact, Confidence, numeric value, tolerance, baseline row deletion, and new success without row. The executable checker must fail every case.
- [ ] Run:

  ```bash
  mamba run -n vntyper pytest -q \
    tests/unit/test_integration_compatibility.py \
    tests/unit/test_integration_outcome_contract.py \
    tests/unit/test_makefile_recipes.py \
    tests/unit/test_workflow_consistency.py
  VNTYPER_TEST_DATA_SKIP_DOWNLOAD=1 mamba run -n vntyper pytest -q \
    tests/integration/test_pipeline_integration.py \
    tests/integration/test_cram_pipeline.py -m integration
  mamba run -n vntyper make check-all
  mamba run -n vntyper make ci-local
  ```

- [ ] Commit:

  ```text
  test(integration): restore and protect real successes

  Closes #233
  ```

---

## Task 6: Gate the fixed behavior in Docker and pin workflow tier selection

**Files:**

- Modify: `tests/docker/test_docker_pipeline.py`
- Modify: `tests/docker/conftest.py`
- Modify: `tests/unit/test_docker_pipeline_runner.py`
- Modify: `tests/unit/test_makefile_recipes.py`
- Modify: `tests/unit/test_workflow_consistency.py`
- Modify: `Makefile`
- Modify: `.github/workflows/docker-build.yml` only if current tier wiring needs correction

### 6A. Add shared-oracle Docker nodes

- [ ] Parametrize Docker BAM, CRAM, and FASTQ from the same declarations and call the same shared validators. Add no Docker-only genotype expectations.
- [ ] Required real nodes: b178 mixed BAM, 40cf 93-singleton BAM, original-derived
  `no_ref=1` mixed CRAM, reference-compressed mixed b178 CRAM, pure-paired remapped CRAM,
  paired b178 FASTQ, and single 6449 FASTQ. Scheduled/full retains mixed a5c1 adVNTR plus
  the pure-paired control.
- [ ] All input/reference mounts are read-only; per-case output starts empty and uses the corrected mapped subdirectory. Assert routing JSON and value/artifact/archive evidence.

### 6B. Make quick a real success sentinel

- [ ] `test-docker-quick` keeps b178 but its declaration must require exit 0, routing singleton, Kestrel/coverage values, summary/report, and archive.
- [ ] Task 6 owns the quick sentinel and records its measured default-mode runtime. It
  must remain below the global 600-second timeout on the two-core Docker runner; if it
  does not, optimize fixture scope rather than moving it off the merge gate.
- [ ] Unit recipe test proves current regression exit 1 cannot satisfy quick.
- [ ] Make `test-docker-fast` and full `test-docker` depend on `cram-fixtures` (quick
  does too if its selected case uses CRAM). Add `docker-cram-fixtures`, which runs
  `scripts/make_cram_fixtures.py` inside the freshly built application image with the
  repository mounted read-only, `tests/data` mounted read-write, host UID/GID, and
  `/opt/vntyper/reference/alignment/chr1.hg19.fa` passed explicitly. This uses the exact
  packaged reference without requiring an absent host `reference/alignment` tree. Fast
  and full Docker targets depend on this target; pin image tag, mounts, entry point,
  reference path, prerequisite, and manifest path in unit recipe tests. The workflow
  invokes it after test-data download and image build but before CRAM parametrization.
- [ ] YAML-semantic workflow contract pins PR→quick, main push→fast, and schedule/manual→full. It asserts the target is executed, not merely present in Make, and no failure masking exists.
- [ ] Any workflow edit requires `make ci-local`; Docker behavior requires `make ci-local-docker`.

### 6C. Build and run exact final-tree Docker matrix

```bash
mamba run -n vntyper make docker-build DOCKER_IMAGE_NAME=vntyper DOCKER_IMAGE_TAG=issue-233
VNTYPER_TEST_DATA_SKIP_DOWNLOAD=1 VNTYPER_TEST_IMAGE=vntyper:issue-233 \
  mamba run -n vntyper pytest -v -m docker \
  'tests/docker/test_docker_pipeline.py::test_docker_bam_pipeline[example_b178_hg19_subset_default]' \
  'tests/docker/test_docker_pipeline.py::test_docker_bam_pipeline[example_40cf_hg38_subset_fast_gdp_guard]' \
  'tests/docker/test_docker_pipeline.py::test_docker_cram_pipeline[b178_hg19_original_derived_mixed]' \
  'tests/docker/test_docker_pipeline.py::test_docker_cram_pipeline[b178_hg19_reference_compressed_mixed]' \
  'tests/docker/test_docker_pipeline.py::test_docker_cram_pipeline[b178_hg19_remapped_pure_paired]' \
  'tests/docker/test_docker_pipeline.py::test_docker_fastq_pipeline[example_b178_hg19_subset_paired_fastq_no_shark]' \
  'tests/docker/test_docker_pipeline.py::test_docker_fastq_pipeline[example_6449_hg19_subset_single_fastq]' \
  'tests/docker/test_docker_pipeline.py::test_docker_advntr_pipeline[example_a5c1_hg19_subset_advntr]' \
  'tests/docker/test_docker_pipeline.py::test_docker_advntr_pipeline[example_b178_hg19_bwa_advntr]'
VNTYPER_TEST_DATA_SKIP_DOWNLOAD=1 VNTYPER_TEST_IMAGE=vntyper:issue-233 \
  mamba run -n vntyper make test-docker-quick
```

- [ ] Mutation proof: drop singleton from Docker routing output, change FASTQ Confidence, skip archive option, or change workflow PR target; each owning test fails.
- [ ] The two explicit Docker adVNTR nodes are mandatory final evidence despite their
  15–25 minute runtime; `ci-local-docker` quick coverage does not substitute for their
  mixed-input Kestrel/adVNTR/cross-match value assertions.
- [ ] Commit:

  ```text
  test(docker): gate mixed-read pipelines end to end

  Closes #233
  ```

---

## Task 7: Correct documentation, verify, review, and publish the PR

**Files:**

- Modify: `docs/superpowers/specs/2026-08-10-milestone-6-harness-matrix-design.md`
- Modify: `docs/plans/2026-08-08-milestone-4-cram-input-robustness-spec.md`
- Modify: `docs/development/golden-cohort-gate.md`
- Modify: `docs/about/changelog.md`
- Modify: `AGENTS.md` if the stable success-baseline trap is adopted
- Add: `tests/unit/test_issue_233_documentation_contract.py`

### 7A. Documentation TDD

- [ ] Add a unit-marked source contract that fails while active docs say `52c41465` introduced refusal, describe valid mixed layouts as intentionally failing, omit no-drop/invalid-parity distinction, or fail to name b178 Docker quick as a success sentinel.
- [ ] Correct history: `2b4597d8` introduced unconditional refusal; `52c41465` inverted nine BAM and one adVNTR oracle. Add dated “superseded by #233” notes rather than rewriting historical evidence.
- [ ] Update golden-cohort documentation from generated fresh evidence only; add changelog wording that valid singleton/unpaired reads are retained without changing clinical text.
- [ ] Run focused contract and strict docs build; commit:

  ```text
  docs(pipeline): record lossless mixed-read policy

  Closes #233
  ```

### 7B. Fresh final verification

From repository root:

```bash
mamba run -n vntyper make format-check
mamba run -n vntyper make lint
mamba run -n vntyper make type-check-all
mamba run -n vntyper make test-unit
mamba run -n vntyper make test-unit-cov
mamba run -n vntyper make test-scripts-cov
mamba run -n vntyper make patch-coverage
mamba run -n vntyper make check-all
mamba run -n vntyper make ci-local
VNTYPER_TEST_DATA_SKIP_DOWNLOAD=1 mamba run -n vntyper make test-integration
VNTYPER_TEST_DATA_SKIP_DOWNLOAD=1 mamba run -n vntyper make ci-local-docker
VNTYPER_TEST_DATA_SKIP_DOWNLOAD=1 mamba run -n vntyper make check-full
```

- [ ] Run the full before/after golden-cohort gate at the final candidate SHA, not only
  its unit tests. Use isolated clean worktrees and retained shared test/reference data:

  ```bash
  python scripts/golden_cohort_gate.py run --side before --tree "$BASE_TREE" \
    --run-root "$RUN_ROOT/before" --marker vntyper.scripts.kestrel_execution \
    --expect-marker absent --threads 2 --jobs 6
  python scripts/golden_cohort_gate.py run --side after --tree "$CANDIDATE_TREE" \
    --run-root "$RUN_ROOT/after" --marker vntyper.scripts.kestrel_execution \
    --expect-marker present --threads 2 --jobs 6
  python scripts/golden_cohort_gate.py compare \
    --before-root "$RUN_ROOT/before" --after-root "$RUN_ROOT/after" \
    --expect-before-sha "$BASE_SHA" --expect-after-sha "$CANDIDATE_SHA" \
    --require-clean --json "$RUN_ROOT/comparison.json" \
    --text "$RUN_ROOT/comparison.txt"
  ```

  Capture every one of the 32 restored mixed cases, all value/artifact deltas, command
  differences, exit codes, and the final verdict. Any genotype delta requires read-level
  attribution; the PR links the retained comparison artifacts.

- [ ] The archive is locally available, so attempt `check-full`. If an external adVNTR dependency fails, record the exact unavailable node and all other passed tiers; do not claim full success.
- [ ] Capture exit codes, counts, branch coverage, scripts coverage, patch coverage, image digest, skips, commands, and artifact paths.
- [ ] Inspect full merge-base diff, marker hygiene, Python 3.10 syntax/runtime, every new failure path, secrets, permissions, and unchanged thresholds.

### 7C. Independent adversarial review

- [ ] Use `superpowers:requesting-code-review` for every task range and whole branch; resolve Critical/Important findings through `receiving-code-review`.
- [ ] Run exact `claude-opus-5`, maximum effort, read-only against issue, design/plan, commits `2b4597d8`/`52c41465`, full diff, real outputs, local/Docker tests, compatibility history, and final logs.
- [ ] Require C/H/M/L and exact evidence. Resolve all Critical/High in one coordinated wave, then one scoped re-review and full affected/final gates.

### 7D. Acceptance traceability

| Requirement | Tasks | Evidence |
| --- | --- | --- |
| AC233-1 lossless mixed routing | 1, 5 | truth table + b178/40cf real runs |
| AC233-2 exactly-once stable inputs | 1, 3, 5, 6 | structured routing + parsed argv + mutations |
| AC233-3 invalid/empty/incomplete fail | 1 | negative unit matrix |
| AC233-4 Kestrel compatibility/quoting | 1, 4 | exact commands + real JAR |
| AC233-5 b178 BAM/CRAM local+Docker | 5, 6 | shared value/artifact oracle |
| AC233-6 downstream input matrix | 3, 5, 6 | FASTQ/CRAM/BAM/adVNTR/cross-match nodes |
| AC233-7 mutation sensitivity | 1–7 | named killed mutations |
| AC233-8 Docker quick success gate | 5, 6 | declaration, Make, workflow, real Docker |
| AC233-9 oracle-inversion protection | 2, 5 | bootstrap, live bidirectionality, event-base tests |
| AC233-10 full gates | 7 | fresh command evidence |
| AC233-11 accurate history/docs | 7 | source contract + strict docs |

### 7E. Push and PR

- [ ] Use `superpowers:verification-before-completion`, `finishing-a-development-branch`, and `github:yeet`.
- [ ] Ensure clean tree, coherent commits, no tracked SDD scratch, no tag operation. Push without force and open PR against `main` with `Closes #233`.
- [ ] PR includes root cause/oracle inversion, architecture, v2.0.3/current/fixed comparison, Kestrel proof, real local/Docker matrix, compatibility guard, mutations, exact gates/coverage, Opus verdict, and risks.
- [ ] Inspect required GitHub CI/Docker checks to terminal green; debug every failure systematically and push scoped follow-ups without force.

## Completion checklist

- [ ] Tasks 1–3 are unit-green and reviewed; Task 1's explicitly bounded stale-oracle
  integration/Docker interval ends at Task 5 and is never presented as a green tier.
- [ ] Atomic production vertical never strands a caller between commits.
- [ ] Equal pairs plus unpaired reads are lossless; unequal/one-sided/empty/incomplete fail.
- [ ] Kestrel 1.0.1 real multi-file behavior is proven.
- [ ] BAM, CRAM, paired/single FASTQ, adVNTR, and cross-match have value-level local/Docker evidence.
- [ ] Docker quick requires mixed b178 success and cannot pass current 2.0.11 behavior.
- [ ] Checker validates live outcomes, bidirectional coverage, historical bootstrap, append-only history, and explicit PR/push bases in CI.
- [ ] All applicable gates and reviews pass, PR exists, and no tag was created or moved.
