# Milestone 6 harness matrix design

**Status:** approved design for the remaining scope of
[#71](https://github.com/hassansaei/VNtyper/issues/71), with
[#226](https://github.com/hassansaei/VNtyper/issues/226) retained as evidence-only
reconciliation.

**Baseline:** `main` at `ebb15b26631242a3295607e4eda4c68f688cd9a2` (`v2.0.10`).

## Purpose

Close the two remaining gaps in #71 without reimplementing the four criteria that the
[criterion-by-criterion audit](https://github.com/hassansaei/VNtyper/issues/71#issuecomment-5214053215)
already found satisfied:

1. pin the positive and negative halves of the intermediate/archive flag matrix with
   objective filesystem assertions; and
2. run a genuinely independent paired FASTQ specimen without SHARK.

The same work corrects the production mismatch between the CLI's documented
`--delete-intermediates` precedence and the implementation. Issue #226 requires no new
fixture work: its missing symbol and reference-resolution coverage were delivered by
milestone 4 and only need to be cited when the issue is closed.

## Source acceptance criteria and disposition

The numbered rows below link to the
[#71 issue acceptance criteria](https://github.com/hassansaei/VNtyper/issues/71).
The body and its audit comment are the source of scope; this design does not broaden
them.

| Criterion | Disposition at the baseline | Exact evidence or required result |
| --- | --- | --- |
| **AC71-1 Docker** | Satisfied; evidence only | `tests/docker/test_docker_pipeline.py` runs parametrised pipelines in the image, and `tests/docker/test_image_structure.py` checks the image contract. No new Docker scenario is required by this design. |
| **AC71-2 API** | Satisfied; evidence only | `tests/unit/web/` covers submission, status, retrieval, validation, authorization and error responses. No API implementation or test change is in scope. |
| **AC71-3 CLI + adVNTR** | Satisfied; evidence only | `integration_tests.advntr_tests` in `tests/test_data_config.json` includes a successful b178 case with `--extra-modules advntr`, `--advntr-max-coverage 300` and value-level checks of `advntr/output_adVNTR_result.tsv`. |
| **AC71-4 Subcommands** | Satisfied; evidence only | Parser, dispatch and handler coverage lives in `tests/unit/test_cli_parser_contract.py`, `test_cli_handlers.py`, `test_cli_dispatch.py`, `test_cli_report.py` and `test_online_mode.py`. End-to-end `report`, `cohort` and `online` runs remain explicitly out of scope. |
| **AC71-5 Parameter matrix** | New work | Successful cases must prove an intermediate is retained with `--keep-intermediates`, removed with `--delete-intermediates`, removed when both flags are present because delete wins, an archive exists with `--archive-results`, and no archive exists when the flag is absent. Assertions must test both presence and absence, not infer behavior from exit code. |
| **AC71-6 Additional FASTQ** | New work | A paired b178 R1/R2 case, distinct from the existing 6449 specimen, must complete without `--extra-modules shark` and must produce the declared success artifacts. The current single-R1 6449 case is useful single-end coverage but is not an independent second specimen. |
| **AC226 fixture symbol** | Satisfied; evidence only | `scripts.make_cram_fixtures.build_reference_dependent_fixture` exists, is called by `main`, and is covered by unit and integration tests. No source change is authorized for #226. |

The fragment identifier in the issue URL is informational rather than a code dependency;
the local test names and paths below are the executable contract.

## In scope

- Declarative presence and absence expectations for successful integration cases in
  `tests/test_data_config.json`.
- Shared output-relative artifact assertions in `tests/support/orchestration.py`.
- Explicit sibling archive presence and absence assertions in
  `tests/integration/test_pipeline_integration.py` and, where the successful matrix is
  hosted, `tests/integration/test_single_end_pipeline.py`.
- Successful derived b178 single-end BAM cases that exercise keep, delete, both flags,
  archive and no-archive behavior without relying on the intentionally failing mixed
  BAM cases.
- A paired `example_b178_hg19_subset_R1.fastq.gz` /
  `example_b178_hg19_subset_R2.fastq.gz` integration case with no SHARK module.
- Correct `--delete-intermediates` precedence in
  `vntyper.scripts.fastq_bam_processing.process_bam_to_fastq`, including a unit test.
- An `Unreleased` changelog entry saying that when both flags are supplied,
  `--delete-intermediates` wins, matching the existing CLI help.
- Unit tests for declaration validation and shared assertion behavior.
- Administrative evidence for closing #226.

## Out of scope

- New Docker, API, adVNTR or subcommand implementations for AC71-1 through AC71-4.
- End-to-end `report`, `cohort` or `online` tests.
- Changing the mixed-layout fail-closed behavior introduced by `52c4146596fef2d1e2402991fbab062ba8021889`.
- Turning the scheduled integration tier into a PR gate.
- Downloading a new fixture: the b178 pair is already declared in
  `tests/test_data_config.json` and the derived single-end BAM already exists in the
  fixture contract.
- Defining every temporary file as public API. Only artifacts named by a test case are
  contractual.
- Any new #226 production code, fixture format or reference-resolution behavior.

## Architecture and boundaries

### 1. Case declarations are the test contract

`tests/test_data_config.json` remains the single source of parametrised integration
inputs. Successful BAM-like cases may declare:

```json
{
  "expected_present": [
    "fastq_bam_processing/output_unmapped.bam"
  ],
  "expected_absent": [
    "fastq_bam_processing/output_unmapped.bam"
  ],
  "expected_archive": true
}
```

The example shows the three fields, not a valid single case: one path must not appear in
both artifact lists.

Interface rules:

- `expected_present` and `expected_absent` are optional lists of POSIX-style paths
  relative to the case output directory.
- Entries must be non-empty, relative, contain no `..` component and be unique.
- The two lists must be disjoint.
- Missing lists mean "no assertion of this kind", not "expect nothing".
- `expected_archive` is optional. When present, its Boolean value is asserted against
  the default sibling ZIP path `Path(f"{output_dir}.zip")`.
- Cases using a non-default archive format must extend the declaration with an explicit
  expected suffix rather than guessing from the filename. No such case is required for
  #71.
- Success-only expectations are forbidden on `expected_exit_code != 0` cases. Extend
  the existing guard in `tests/unit/test_integration_outcome_contract.py` to include
  `expected_present` and `expected_absent`.

The paths are deliberately explicit. A broad assertion such as "the output directory is
smaller" would not say which behavior changed, and a directory-exists assertion could
pass while the load-bearing intermediate remained behind.

### 2. Shared orchestration owns output-relative assertions

Add one focused helper in `tests/support/orchestration.py`, conceptually:

```python
assert_declared_artifacts(test_case: dict, output_dir: Path) -> None
```

It validates the declaration, resolves every entry under `output_dir`, and reports all
missing and unexpectedly present paths in one assertion. `run_bam_test_case` calls it
only after `_assert_expected_exit` confirms a successful application outcome and after
the ordinary required-file checks. Unit tests drive it with `tmp_path`; it performs no
network or subprocess work.

Sibling archives do not live under `output_dir`, so their assertion stays in the local
integration caller. Replace truthiness-only logic with key-presence logic:

- `expected_archive: true` requires the sibling ZIP;
- `expected_archive: false` requires that the sibling ZIP does not exist;
- an omitted key makes no archive assertion.

This keeps `run_bam_test_case` usable by the current Docker runner, which does not pass
the case's arbitrary CLI options. AC71-5 is exercised by the local scheduled integration
tier, not smuggled into Docker as a different command.

### 3. Use a successful homogeneous BAM path for the matrix

The nine current `integration_tests.bam_tests` stop at the mixed-read-layout refusal and
all expect exit 1. They cannot reach cleanup or archive creation and must not be
repurposed as false success cases.

Host the matrix on `integration_tests.single_end_bam_tests`, whose b178 fixture is
homogeneous and currently completes. The minimum declarations are:

| Case | CLI options | Required assertion |
| --- | --- | --- |
| keep | `--keep-intermediates`, no `--fast-mode` | `fastq_bam_processing/output_unmapped.bam` present; sibling ZIP absent |
| explicit delete | `--delete-intermediates`, no `--fast-mode` | `fastq_bam_processing/output_unmapped.bam` absent |
| precedence | `--keep-intermediates --delete-intermediates`, no `--fast-mode` | `fastq_bam_processing/output_unmapped.bam` absent |
| archive | `--archive-results` | sibling ZIP present after exit 0 |

Cases may be combined only when every row remains independently diagnosable. In
particular, the precedence case must carry both flags and its own ID; a delete-only case
cannot prove the documented conflict rule.

`--fast-mode` is intentionally absent from the cleanup rows because that mode never
creates the unmapped BAM whose retention is being tested.

### 4. Production precedence is enforced at the cleanup owner

The parser already says `--delete-intermediates` "overrides
--keep-intermediates" in `vntyper/scripts/cli_parser.py`. The current cleanup condition
inside `process_bam_to_fastq` is the inverse for the conflict case:

```python
if delete_intermediates and not keep_intermediates:
```

`pipeline.py` always passes `output_name="output"` to `process_bam_to_fastq`, so the
literal artifact under test is `fastq_bam_processing/output_unmapped.bam`; fixture or
case names must never be interpolated into that assertion. The keep-row RED test first
proves this path is real, and only then may the two absence rows use the same literal.
This sequencing prevents a negative assertion against an artifact the pipeline never
creates from satisfying #71.

Change the cleanup owner's decision so an explicit delete request removes the declared
intermediate even when keep is also true. Update the function's argument documentation
to state the precedence. A direct unit test must fail against the old condition and
prove all four Boolean combinations, with particular emphasis on `(delete=True,
keep=True) -> delete`.

The behavior change is documented above the released 2.0.10 section in a new
`Unreleased` changelog section. Released history is not rewritten.

### 5. The second FASTQ must be independent

Add a paired case using:

- `tests/data/example_b178_hg19_subset_R1.fastq.gz`
- `tests/data/example_b178_hg19_subset_R2.fastq.gz`
- reference assembly `hg19`
- no `--extra-modules` entry
- the normal success artifacts, including `summary_report.html` and
  `kestrel/kestrel_result.tsv`.

The case must assert both input paths and that SHARK is absent from its parsed module
list in a pure contract test. Reusing 6449 R1 as a single-end case remains useful for
single-input routing but does not count as AC71-6.

### 6. #226 is an evidence boundary, not an implementation task

No implementation file changes for #226 are part of this track. Closure evidence is:

- `build_reference_dependent_fixture` at `scripts/make_cram_fixtures.py:340`, called by
  `main` at line 464;
- introduction in `80f42fbba63c5c045ee605e49949cdd970884320` and hardening in
  `edabd3eaf594b785906d8ba03c9cce60f1c6babd`;
- `tests/unit/test_make_cram_fixtures.py::test_reference_dependent_fixture_has_a_local_ur_target_that_can_be_removed`;
- `tests/integration/test_cram_reference_contract.py::test_a209_1_missing_reference_names_the_digest_and_candidates_before_stages`;
- `test_a209_2_explicit_reference_completes_the_reference_dependent_cram`; and
- `test_a209_3_no_ref_cram_completes_without_an_explicit_reference`.

The issue is closed administratively with those links after the harness PR; it receives
no `Closes #226` code change manufactured after the fact.

## Data flow and invariants

1. Pytest loads the case from `tests/test_data_config.json` at repository-root CWD.
2. The integration runner creates a fresh case-specific `tmp_path` output directory.
3. The runner appends the case's CLI options without reordering them.
4. VNtyper exits. Non-zero expected outcomes stop before success-only assertions.
5. On exit 0, shared orchestration validates required reports, values and declared
   output-relative artifacts.
6. The local integration caller validates declared sibling archive state.
7. Every case gets a fresh output directory, so an artifact from another case cannot
   satisfy an assertion.

Hard invariants:

- `--delete-intermediates` wins whenever both keep and delete are true.
- An absence check is a first-class assertion, never inferred from a missing positive
  declaration.
- Artifact paths cannot escape the case output directory.
- Failure cases cannot promise unreachable success artifacts.
- The b178 FASTQ case is paired, uses a specimen other than 6449 and omits SHARK.
- #226's existing fixture and reference tests are not altered to make #71 pass.

## Success, partial failure, recovery, retry and idempotency

### Success

A case succeeds only when the application exit matches, ordinary result validation
passes, every `expected_present` path exists, every `expected_absent` path does not
exist, and the sibling archive state matches when declared.

### Partial failure

- A pipeline failure reports the command output through the existing integration
  runner and does not run success-only artifact assertions.
- A declaration error fails before filesystem assertions with the case ID and invalid
  field/path.
- Artifact mismatches list the exact resolved paths and whether each was missing or
  unexpectedly present.
- Archive creation failure is not accepted merely because the output directory itself
  is valid.

### Recovery and retry

All mutable test output lives under a fresh pytest temporary directory. A failed or
interrupted case does not alter test data and requires no recovery command. Re-running
the node creates a new output directory and evaluates the same declarations.

### Idempotency

The test contract is idempotent because inputs are read-only, output roots are unique,
and declarations contain no timestamps or generated paths. Production archive safety
continues to own overwrite/refusal behavior; this track does not weaken it for tests.

## Compatibility and environment

- Implementation and tests must run on Python 3.10 through 3.13.
- New unit test files carry `pytestmark = pytest.mark.unit` and use only `tmp_path` and
  mocks.
- Integration pytest is launched from the repository root because
  `tests/parametrization.py` reads `tests/test_data_config.json` by relative path.
- The b178 inputs and derived fixture come from the existing test-data contract; no
  network is added to unit tests.
- Integration data bootstrap remains fail-closed. Use
  `VNTYPER_TEST_DATA_SKIP_DOWNLOAD=1` when a missing archive should fail immediately.
- The integration tier remains scheduled rather than a PR gate. That limitation must
  be stated in the PR and issue closure evidence.

## Security and permissions

- Test-declared relative paths are confined under `output_dir`; absolute and parent
  traversal paths are rejected before access.
- Inputs are never deleted or rewritten. Cleanup assertions name only files created
  below the case output root.
- No credentials, GitHub permissions, network endpoints or workflow token scopes change.
- No Docker container is granted additional mounts or write access for this track.

## Observability and diagnostic output

- Pytest IDs retain the case `test_name`, including explicit `keep`, `delete`,
  `delete_overrides_keep`, `archive`, and `b178_no_shark` wording.
- Failure messages name the case, declaration field and resolved artifact.
- Integration logs continue to capture VNtyper stdout and stderr.
- The precedence changelog entry is the user-visible record of the corrected behavior.

## Test strategy and objective acceptance checks

| Test | Proves |
| --- | --- |
| Unit test of declaration validation | Paths are relative/confined, lists are disjoint, and invalid declarations fail closed. |
| Unit test of `assert_declared_artifacts` | Both missing-required and unexpectedly-present artifacts fail with exact paths. |
| Unit test of negative-case hygiene | Exit-1 cases cannot retain `expected_present`, `expected_absent` or `expected_archive`. |
| Direct unit test of cleanup precedence | All flag combinations are pinned and both flags select deletion. |
| Successful single-end keep integration node | AC71-5 positive intermediate half. |
| Successful single-end explicit-delete node | AC71-5 negative intermediate half. |
| Successful single-end both-flags node | Documented delete-overrides-keep contract. |
| Successful archive/no-archive nodes | AC71-5 positive and negative sibling archive halves. |
| Pure FASTQ declaration contract test | b178 R1 and R2 are used and no SHARK option is present. |
| Paired b178 integration node | AC71-6 end-to-end alternate specimen without SHARK. |
| Existing #226 unit and A209 integration nodes | #226 evidence-only acceptance; no new implementation. |

Required verification:

```bash
make test-unit
make test-unit-cov
make patch-coverage
make test-integration
make check-all
```

If Docker-facing orchestration is changed beyond the boundaries above, additionally run
`make ci-local-docker`; it is not required merely to prove the already-satisfied
AC71-1.

## Rollout

1. Add pure declaration/assertion tests first and show they fail against the absent
   schema behavior.
2. Correct production precedence with its direct regression test and add the
   `Unreleased` changelog entry.
3. Add the successful single-end matrix declarations and shared assertions.
4. Add the independent paired b178 no-SHARK declaration and contract test.
5. Run unit gates, then the data-dependent integration tier.
6. Close #71 with test-node and command evidence for AC71-5 and AC71-6 only.
7. Close #226 separately with the existing symbol, commits and A209 tests cited above.

## Unresolved questions

None. The remaining issue scope, delete-overrides-keep precedence, artifact schema,
successful fixture host and independent FASTQ specimen are fixed by this design.
