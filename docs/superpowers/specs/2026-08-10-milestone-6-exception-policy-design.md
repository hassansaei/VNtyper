# Milestone 6 exception-policy design

**Status:** Approved implementation specification

**Date:** 2026-08-10

**Baseline:** `main` at `ebb15b26631242a3295607e4eda4c68f688cd9a2`

**Issue:** [#219](https://github.com/hassansaei/VNtyper/issues/219)

**Milestone:** 5. Gates, harness and release automation (parallel)

## Decision summary

Keep `BLE001` out of Ruff's global selected rules. Replace the stale prose-only exemption with
two exact, executable aggregate guards and a symbol-based audit of handlers that swallow an
exception and continue with a fallback, empty, negative, unavailable, or partial result.

At the specified baseline, clean no-cache measurements are:

- **103** BLE001 diagnostics under normal Ruff suppression semantics;
- **108** BLE001 diagnostics with `--ignore-noqa`, exposing the five existing explicit
  suppressions.

The implementation must remeasure both numbers at its final tip and pin the exact final pair.
The expected values begin at 103/108; they change only when a reviewed handler change or explicit
suppression in this PR explains the delta. The old issue-era counts of 46 and 67 are historical
and must not appear as the current inventory.

Every fail-open handler is identified by module and qualified symbol rather than a source line,
assigned a policy disposition, given a concise domain rationale, and linked to a behavior test.
An intentional fail-open result is preserved only when that rationale and test exist. Behavior
changes are allowed only when an existing contract supplies the replacement outcome. Where no
alternate outcome is authorized, the present behavior is characterized and preserved rather
than guessed.

This incorporates the two review corrections: aggregate growth must be caught exactly under both
suppression modes, while the high-risk audit must be stable across line movement and focus on
semantic fail-open behavior rather than mechanically rewriting every broad catch.

## Linked issue and exact acceptance criteria

Issue #219 says the current all-or-nothing non-selection is reasonable and explicitly does not
request reversing it wholesale. Its acceptance criteria, made objective here, are:

1. Correct the stale BLE001 count and make it self-checking so growth is visible.
2. Classify the inventory into:
   - genuine process or stage boundaries where catch-all plus a controlled failure outcome is
     the contract;
   - handlers around calls whose exception surface is small enough to enumerate, but which need
     call-specific evidence before narrowing;
   - fail-open handlers that swallow an exception and continue with a degraded or neutral result.
3. Give fail-open handlers priority because they can make an operational failure look like a
   valid clinical result.
4. Do not turn BLE001 on globally.
5. Do not narrow the complete inventory in one pass; each narrowing requires evidence about the
   called API and the enclosing contract.
6. Leave `G004` and the repository's established logging-string style out of scope.

This specification adds the following review-mandated acceptance criteria:

7. Pin both the normal count and the `--ignore-noqa` count, not a raw grep count and not only one
   Ruff view.
8. Map every counted diagnostic to exactly one classification using source/AST information.
9. Maintain a separate, stable-symbol manifest for every fail-open classification. Line numbers
   are observations for error messages, never identities.
10. Each preserved fail-open behavior has a rationale and an objective behavior test.
11. A behavior change uses an already documented or tested failure outcome; the audit does not
    invent a new exit code, fallback, exception class, clinical message, or pipeline policy.

## Scope

### In scope

- All BLE001 diagnostics emitted from the same path set as `RUFF_PATHS`:
  `vntyper/`, `docker/app/`, `tests/`, `scripts/`, and `docs/`.
- Both Ruff suppression modes: normal operation and `--ignore-noqa`.
- `except Exception` and `except BaseException` handlers that Ruff reports as BLE001.
- A complete inventory classified by enclosing module/qualified symbol and policy category.
- A focused fail-open audit covering every category-C symbol, including at minimum the current
  high-signal production symbols:
  - `vntyper.scripts.flagging.regex_match`;
  - `vntyper.scripts.flagging.evaluate_condition`;
  - `vntyper.scripts.variant_parsing.read_vcf_without_comments`;
  - `vntyper.scripts.pipeline_guards.read_alignment_header`;
  - `vntyper.scripts.region_utils._bam_contig_names`;
  - `vntyper.scripts.report_formatting.parse_coverage_stats`;
  - `vntyper.scripts.screening_summary.load_report_config`;
  - `vntyper.scripts.screening_summary.build_screening_summary`;
  - `scripts.coverage_gate.read_total`.
- Any additional fail-open handler found by the complete classification. The minimum list is not
  an allowlist and does not cap the audit.
- Tests of documented fallback, fail-closed, exit-status, re-raise, and unavailable-result
  contracts.
- Updating the BLE001 rationale comment with the final measured pair and a pointer to the policy
  test, while retaining explicit Ruff `select` semantics.
- Contributor documentation and the MkDocs exclusion of `docs/superpowers/` working
  specifications, owned by the quality-gates rollout.

### Out of scope

- Adding BLE001 to `[tool.ruff.lint].select` globally.
- Adding blanket `per-file-ignores`, mass `# noqa: BLE001` comments, or a repository-wide baseline
  file that silently accepts future growth.
- Mechanically changing `except Exception` to a guessed tuple of exception classes.
- Changing a return, exit, retry, or fallback contract solely because a linter diagnostic exists.
- Reworking the already fail-closed `cross_match.match_variants`; commit `5ee1e4a` is a design
  precedent, not a target for reversal.
- Changing flag names, stage filtering, report wording, clinical classification, or pipeline exit
  codes without an existing source of truth that specifies the replacement outcome.
- G004, logging format conversion, custom exception classes, or removal of intentional stage and
  process boundaries.
- Network, Docker, external tools, or reference data in the unit-tier policy tests.

## Considered approaches

### Approach A: select BLE001 globally and suppress boundary files

This is rejected. It would turn a semantic inventory into a large suppression exercise, obscure
which catches are safe boundaries, and make file-level suppressions cover future handlers. It
also conflicts with the repository's explicit-rule policy and with the issue's stated intent.

### Approach B: narrow every broad catch to obvious built-in exceptions

This is rejected. External binaries, libraries, parsers, templates, filesystem operations, and
process boundaries have different exception surfaces. Guessing a tuple can turn a controlled
failure into an unattended crash, while a superficially narrower catch says nothing about whether
the fallback is safe.

### Approach C: pin counts, classify all diagnostics, and test fail-open contracts by symbol

This is selected. Exact dual counts stop silent aggregate growth, full classification records why
the current non-selection exists, and symbol-based fail-open tests protect the clinically relevant
semantics across routine refactors. Narrowing remains possible one handler at a time when evidence
supports it, without making narrowing itself a milestone success metric.

## Policy taxonomy and decision rule

Every diagnostic belongs to exactly one category:

### A. Controlled boundary

The handler is at a CLI, task, stage, process, report-generation, or external-command boundary and
converts arbitrary failure into a defined terminal outcome: re-raise, `parser.error`, `sys.exit(1)`,
a failed task state, an HTTP error, or another explicit non-success result. The audit records the
terminal outcome. Catching broadly remains permitted when the boundary genuinely owns all failure
translation.

### B. Enumerable candidate

The handler surrounds a local or documented third-party call whose exception set may be narrowed,
but the current milestone has not yet demonstrated the complete set. The handler remains unchanged
unless tests, API documentation, and caller behavior jointly prove the replacement catch list.
Classification as B is not authorization for a mechanical edit.

### C. Fail open or degraded continuation

The handler swallows an exception and returns or continues with a value that downstream code can
mistake for ordinary output: `False`, `None`, an empty collection/frame/config, a negative call,
an unavailable summary, default metrics, skipped validation, or a partial artifact. Every C entry
must be present in the stable-symbol audit with a disposition, rationale, and behavior test.

Use this deterministic disposition rule for each C entry:

1. If public documentation, configuration semantics, an established test, or an explicit caller
   contract names the error outcome, use that exact outcome. If current code differs, a focused
   behavior change may make it conform.
2. If the existing contract explicitly defines the fallback as meaningful degradation, preserve
   it and add a rationale plus a behavior test that distinguishes failure from success in logs or
   structured state.
3. If no source of truth defines an alternative, preserve the observed behavior for this
   milestone, characterize it, and record the disposition as `preserved-no-authorized-alternative`.
   This is a final policy decision, not a deferred design question: the audit gains visibility
   without inventing domain behavior.

The audit never treats “the function did not raise” as sufficient. A test asserts the exact
return/exit/state, log severity and message identity where stable, and whether downstream control
flow continues or stops.

## Architecture and component boundaries

### Ruff measurement adapter

A unit-test helper invokes the repository-installed `ruff` executable twice with `--no-cache`,
`--select BLE001`, and the exact `RUFF_PATHS` scope:

1. normal suppression semantics;
2. the same command with `--ignore-noqa`.

Machine-readable JSON output is parsed for inventory and human-readable statistics are retained
as the documented reproduction command. The adapter accepts the repository root and Ruff path as
inputs and returns normalized diagnostics containing relative path, row, column, code, and message.
It never edits source or executes imported application functions.

The checked policy metadata records the Ruff version used for the reviewed baseline
(`ruff 0.14.3` at reconstruction). The adapter always records the actual `ruff --version`
alongside its diagnostics. A version change alone is not a failure when the normalized
diagnostic inventory is identical. If identities or counts differ, the assertion names
both reviewed and actual versions and requires an intentional inventory reclassification;
it never reports a source-growth regression without attributing the measuring tool. This
keeps the exact no-unclassified-handler contract while avoiding a blind raw-count failure.

### Classification inventory

An AST mapper assigns each Ruff diagnostic to its enclosing module and qualified function or
class-method symbol. Multiple diagnostics in one symbol are represented by an expected count for
that symbol; module-level handlers use the explicit `<module>` symbol. The inventory records:

- relative source path;
- qualified symbol;
- expected number of BLE001 diagnostics;
- category A, B, or C;
- one-sentence classification rationale.

The inventory covers the 108 `--ignore-noqa` diagnostics, so explicit suppressions are classified
rather than invisible. The normal inventory must be an exact subset and total the separately
pinned normal count. A diagnostic that cannot be mapped, maps to two entries, or is absent from
the inventory fails the test.

The reviewed schema separately pins 79 stable identities and identity-category totals A=30, B=16,
and C=33. These are checked independently of diagnostic totals, so a category-only A-to-B change
cannot pass merely because the same handlers and 103/108 diagnostic counts remain.

Rows and columns appear only in failure output. Renaming or moving a symbol requires an explicit
inventory change; merely adding blank lines does not.

### Stable-symbol fail-open audit

Category C is additionally represented by a manifest keyed by
`<module path>:<qualified symbol>`. Each entry contains:

- disposition: `preserved-contract`, `conformed-to-existing-contract`, or
  `preserved-no-authorized-alternative`;
- fallback or terminal outcome in structured terms;
- domain rationale;
- behavior-test node ID;
- whether the failure is externally observable through a log, status, exception, or result field.

A policy test resolves every symbol from source/AST, confirms that it still contains the classified
broad handler, verifies that each C inventory symbol has exactly one audit entry, and rejects audit
entries for symbols no longer classified C. This protects against line drift and against a handler
quietly moving to an unaudited helper.

The manifest is hand-reviewed policy data, not generated from current code. Ruff/AST discovery is
generated; rationale and disposition are not.

Behavior-test IDs may name top-level or class-qualified tests under any nested `tests/unit/`
subdirectory, including parametrized suffixes. Paths outside `tests/unit/`, parent traversal,
missing `test_` function components, and newline-bearing IDs are invalid. The source-read-only live
validator resolves every frozen node immediately: before Tasks 4–9 create those nodes it exits 1
and enumerates every unresolved ID, and it cannot report complete success until all IDs resolve.

### Behavior tests

Behavior tests live beside the affected module's existing unit tests and use the current unit
marker. They drive the symbol with a deterministic exception from a mocked dependency, then assert
the contractual outcome. Existing tests are extended rather than duplicated where they already
characterize the behavior.

Each node is first run unchanged as characterization. After it passes, a temporary behavior mutation
must make the exact value/state assertion fail and must be restored before policy validation. A
passing characterization is not mislabeled as an expected RED implementation failure.

Important seed contracts include:

- `flagging.evaluate_condition`: test missing names and other evaluation errors separately; assert
  the current `False` gate result and warning/error visibility unless an existing configuration
  contract explicitly requires a different outcome.
- `variant_parsing.read_vcf_without_comments`: distinguish a valid empty VCF, a missing file, and
  an unexpected read/parse failure even if the present DataFrame fallback remains the same.
- `pipeline_guards.read_alignment_header` and `region_utils._bam_contig_names`: assert their
  documented `None`/empty evidence outcomes and the caller's subsequent assembly/region decision.
- `report_formatting.parse_coverage_stats`: assert the complete default metric structure and error
  log rather than only an empty/non-empty result.
- `screening_summary` loaders/builders: assert that failures produce the established unavailable
  state and never a positive emphasis inferred from message text.
- `generate_report.extract_igv_content`: an unreadable report returns the exact three-value tuple
  `("", "", "")`, matching its unpacking caller; it is not a scalar empty-string fallback.
- `coverage_gate.read_total`: assert `None` is converted by the caller into a non-zero gate result
  with an actionable message, so this is degraded parsing but not a fail-open CI pass.

### Configuration ownership

`pyproject.toml` remains the authority for Ruff rule selection and explanatory comments. The policy
test is the authority for exact counts and classification. The Makefile's `RUFF_PATHS` remains the
scope authority. Its tokens must be existing ordinary repository-relative paths, and every Ruff
command places `--` before them so an option-shaped token cannot become an output or mutation
option. A test confirms that policy measurement uses the same paths rather than a duplicated list.

## Interfaces, inputs, outputs, and invariants

| Interface | Inputs | Outputs | Invariants |
| --- | --- | --- | --- |
| normal BLE001 measurement | repository root, Ruff executable, `RUFF_PATHS` | Ruff version, normalized diagnostics and exact count | Baseline is 103; final inventory equals a clean remeasurement. A mismatch reports reviewed and actual Ruff versions. |
| all BLE001 measurement | same plus `--ignore-noqa` | Ruff version and normalized diagnostics including suppressed handlers | Baseline is 108; count is never below normal; the delta equals explicitly suppressed diagnostics. |
| AST classifier | diagnostic path/row and parsed source | module, qualified symbol, category entry | Every all-mode diagnostic maps once; no line number is a persistent identity. |
| fail-open manifest | category-C symbols | disposition, outcome, rationale, test node ID, observability | Exact one-to-one coverage of C; no empty rationale or missing behavior test. |
| behavior test | deterministic dependency exception | asserted fallback/terminal result plus visible signal | Test fails if failure becomes ordinary success, silently changes outcome, or bypasses logging/state. |
| Ruff rationale comment | final measured counts and policy-test path | contributor guidance | States that BLE001 is deliberately not selected; does not claim every handler is a boundary. |

Global invariants:

- `BLE001` remains absent from the global Ruff `select` list.
- No broad per-file ignore is introduced.
- Normal and `--ignore-noqa` counts are exact equality guards, not upper bounds.
- The complete classified total equals the all-mode count; A + B + C equals that same total.
- Every C symbol has a reviewed rationale and behavior test.
- A count update without a corresponding inventory and rationale change fails review and tests.
- Existing completed-application exit codes remain 0/1 and argparse usage errors remain 2.
- Clinical-sounding output text is not invented or reworded.

## Dependency and implementation ordering

### Phase 1: reproduce and pin the dual baseline

1. From a clean tree, run both no-cache Ruff commands over the exact Makefile scope.
2. Confirm the baseline pair is 103 normal and 108 with `--ignore-noqa`.
3. Add exact aggregate tests and a scope-parity assertion while leaving production handlers
   unchanged.
4. Update the pyproject rationale from the stale 46/67-era language to the remeasured pair and
   point contributors to the executable inventory.

### Phase 2: classify the complete all-mode inventory

1. Map every one of the 108 baseline diagnostics to a stable enclosing symbol.
2. Apply categories A, B, and C using the definitions above.
3. Assert that category totals equal the all-mode count and that the normal set is the exact
   unsuppressed subset.
4. Fail the phase if any diagnostic is unclassified, multiply classified, or hidden by a path not
   present in `RUFF_PATHS`.

### Phase 3: audit fail-open behavior

1. Create one complete manifest record for every category-C symbol found by Phase 2, including
   disposition, exact outcome, rationale, behavior-test node ID, and observability.
2. Run the full live validator while linked behavior nodes are being added; require exit 1 with an
   exhaustive unresolved-node list, and do not claim successful completion while any node is absent.
3. Add or strengthen every linked behavior test, run unchanged characterization, and prove its exact
   assertion fails under a temporary mutation before restoring the source.
4. Require the already-enabled full live validator to pass only after every frozen C record's
   behavior node resolves.
5. For an existing contract that requires a different outcome, make only the focused semantic
   change needed to conform and test both the failure outcome and its caller-visible effect.
6. Preserve every other current fallback explicitly; do not narrow category B as a side quest.

### Phase 4: final remeasurement and integration

1. Re-run both no-cache Ruff measurements after all reviewed semantic edits.
2. Set the exact expected pair to the final results. Any delta from 103/108 must be attributable to
   named inventory changes in this PR; an offsetting addition and removal still require two
   explicit classification changes.
3. Re-run the classification, symbol, behavior, lint, type, unit, coverage, and docs gates.
4. Ensure `docs/superpowers/` is excluded from MkDocs through the quality-gates rollout because it
   is contributor policy material, not published user documentation.

## Success, partial failure, recovery, retry, and idempotency

### Success

The work succeeds when both Ruff modes equal their final pinned counts; all all-mode diagnostics
are classified exactly once; every category-C symbol has a disposition, rationale, observable
outcome, and passing behavior test; global BLE001 selection remains off; and repository gates pass.
Reducing the count is not required for success. Making risk visible and preventing silent growth is
the milestone outcome.

### Partial failure

- If the normal count differs but the all-mode count does not, inspect added or removed `noqa`
  suppression. Do not update a number until the suppression is explicitly classified.
- If both counts differ, print diagnostic additions/removals by relative path and symbol. Update
  inventory and counts only for reviewed source changes; otherwise fix the accidental drift.
- If a diagnostic has no stable enclosing symbol, classify it under `<module>` rather than using a
  line number.
- If a C symbol lacks an existing replacement outcome, preserve and characterize its current
  behavior under `preserved-no-authorized-alternative`; do not invent a failure mode.
- If a proposed narrowing causes an unhandled exception or changes a caller outcome, revert that
  narrowing while retaining the classification and characterization test.
- If a behavior test requires network, Docker, reference data, or a real subprocess, extract or
  mock the effect boundary; it does not move to the integration tier.

### Recovery and retry

The Ruff adapter is read-only and uses `--no-cache`, so repeated measurements are independent of
prior runs. AST and manifest checks read tracked source only. Behavior tests use fresh fixtures and
restore patched globals/environment through pytest fixtures. After a failure, correct the source,
inventory, or test indicated by the structured diff and rerun the failed phase followed by both
aggregate measurements.

### Idempotency

Repeated audit execution does not modify source, Ruff configuration, caches, or external state.
Normalized diagnostic ordering is deterministic by relative path, symbol, row, and column for
reporting, while equality is computed from normalized records rather than terminal ordering.
Temporary files exist only below pytest-managed directories.

## Compatibility and environment constraints

- The audit and all behavior tests run on Python 3.10 through 3.13. AST handling and annotations
  use only Python 3.10 APIs and syntax.
- The Ruff version installed by the development/CI environment is the measuring implementation;
  tests always pass `--no-cache` and the explicit rule code.
- Tests run from the repository root with `pytest -m unit`, and each new test module declares
  `pytestmark = pytest.mark.unit`.
- The policy test derives path scope from the Makefile and compares it to the expected repository
  areas. It does not depend on shell glob expansion or the caller's current directory.
- No test requires bwa, samtools, fastp, bcftools, Java, adVNTR, Docker, or the Zenodo archive.
- If workflow files change, the explicit uv virtual-environment install path is verified through
  `make ci-local`.

## Security and permissions

- Ruff and AST inventory operations are source reads only. They do not import, eval, or execute
  scanned repository code.
- Tests of `evaluate_condition` use fixed synthetic expressions and rows. They do not broaden its
  eval namespace or introduce untrusted input into test helpers.
- Mocked exception messages contain no sample data, tokens, paths outside temporary fixtures, or
  environment secrets.
- The change requires no additional GitHub Actions permissions, secrets, registry access, network
  access, or external service calls.
- Failure logs report relative repository paths and stable symbols, not sensitive absolute paths
  or genomic input content.

## Observability and diagnostic output

An aggregate mismatch reports:

- command mode (`normal` or `ignore-noqa`);
- expected and actual count;
- added and removed normalized diagnostics grouped by stable symbol;
- suppression delta and the explicit-suppression symbols responsible for it;
- category totals A/B/C.

A fail-open audit mismatch reports the missing or stale module-qualified symbol, its current
category, disposition, linked behavior-test node ID, and absent field. Behavior tests assert the
existing logger severity or structured error state so preserved fallbacks are observable rather
than silent.

The pyproject rationale reports the final exact pair and explains that the larger number is the
complete inventory including explicit suppressions. It does not use raw grep totals, which count
handlers Ruff intentionally excludes.

## Test strategy and objective acceptance checks

The following checks are mandatory and conjunctive:

1. A clean normal command over `vntyper/ docker/app/ tests/ scripts/ docs/` reproduces the pinned
   final normal count (103 before reviewed implementation deltas).
2. The same command with `--ignore-noqa` reproduces the pinned final all-mode count (108 before
   reviewed implementation deltas).
3. A mutation that adds one unsuppressed blind catch fails both aggregate and classification
   tests; adding one suppressed catch fails the all-mode and classification tests.
4. A mutation that removes or renames a category-C symbol fails the stable-symbol audit without
   depending on its line number.
5. A mutation that removes a C rationale, disposition, observable outcome, or behavior-test link
   fails policy validation.
6. A behavior mutation that turns a dependency exception into ordinary success fails the linked
   module test.
7. Inventory assertions prove all-mode diagnostics equal A + B + C and every diagnostic maps
   exactly once.
8. Configuration assertions prove BLE001 is absent from global `select`, no broad new per-file
   ignore exists, and the measurement scope equals `RUFF_PATHS`.
9. `make lint`, `make type-check-all`, `make test-unit`, and `make test-unit-cov` exit 0.
10. `make check-all` exits 0; `make ci-local` also exits 0 if workflows change.
11. `make docs-build` exits 0 and does not publish or template files below
    `docs/superpowers/`.

## Rollout

Land the policy in one PR using reviewable, green commits:

1. exact 103/108 baseline guards, Ruff-scope parity, and corrected pyproject rationale;
2. complete A/B/C inventory keyed by stable enclosing symbols;
3. category-C manifest plus characterization tests for preserved behavior;
4. any individually authorized conformance changes, each paired with its caller-visible behavior
   test;
5. final dual remeasurement, exact constants, contributor documentation, and full verification.

The change has no runtime feature flag and no data migration. Most of the rollout is test and
policy enforcement. A focused runtime behavior change is independently revertible without
discarding the inventory or tests for other symbols. After merge, a contributor adding or moving a
broad catch must update the exact count, classification, and—if it fails open—the symbol audit and
behavior test in the same PR.

This specification contains no deferred design decisions. The rule selection, count semantics,
classification algorithm, fail-open disposition rule, and allowed behavior-change boundary are
fixed above.

## Unresolved questions

None.
