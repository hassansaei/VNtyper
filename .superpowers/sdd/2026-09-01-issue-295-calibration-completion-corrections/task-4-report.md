# Task 4 report: fail-closed custody and CLI outcomes

## Status

DONE_WITH_DOCUMENTED_LIMITATION

Implementation commit:

- `a3731e7030912332811aff4eae8c61039c8da699 fix(calibration): enforce external custody bindings`

No Task 5, golden, release, version, tag, workflow, or issue-administration work was performed.

## Outcome

Ordinary extraction no longer writes locked-heldout features, labels, baseline values,
plaintext locked payload, or the former locally generated external-custodian evidence
manifest. The locked role contains only a value-free member declaration, exact run-root
and artifact-hash commitments, and its direct-child checksum manifest.

`evaluate --evidence` now treats its input as a separately supplied, closed four-file
custodian import. It rejects an ordinary local extraction bundle and requires:

- a strict named external-custodian authority assertion;
- a passed validation attestation bound to the explicit profile and protocol;
- exact study, protocol, partition, fitted-profile dataset, locked dataset, locked
  run-artifact, locked payload-byte, and validation-attestation-byte hashes; and
- exact inventory, regular non-symlink files, and closed checksum fields.

The local verifier deliberately does not claim that those checks prove independence.
There is no configured cryptographic trust root in this task. The output records that it
verified a supplied named authority assertion and exact bindings, while explicitly
stating that this is not cryptographic proof of external independence. Production has
no command that mints a custodian import; the fixture builder exists only in unit tests.

## Changed contracts

- Added a closed study binding for generated profiles. Its dataset commitment includes
  the complete study, protocol, partition, objective, seed, training/policy-selection
  role manifests, and their exact run-artifact commitments.
- Validation checks the generated profile against that binding before evaluation.
- Held-out evaluation requires the external authority to bind the same profile/study
  lineage plus exact locked run and payload commitments.
- Locked evidence receives a durable precommit before payload bytes are opened and a
  one-use consumption receipt before read access. Exceptional open/evaluation paths
  install a deterministic retirement record and re-raise.
- Completed failed validation and held-out evaluation write a failed attestation and
  retirement record and return `False`. Atomic CLI installation happens before `False`
  becomes exit status 1. Completed success returns `True` and keeps status 0; argparse
  usage remains status 2.
- Retired profile/evidence pairs are rejected before another operation. Retirement is
  pair-specific; it does not retire unrelated profiles.
- Atomic staging and exclusive custody writes catch `BaseException` only at their exact
  cleanup/retirement boundaries, remove only their own partial target, and re-raise
  `KeyboardInterrupt` or `SystemExit`.
- Extracted the new pure locked-artifact, result-encoding, study-binding, and custodian-
  authority logic. `calibration_artifacts.py` is 635 lines after the extraction, below
  the repository's approximate 650-line guideline.

The command surface remains exactly `extract`, `fit`, `validate`, and `evaluate`. The
packaged/default neutral profile was not changed.

## RED-first evidence

The initial semantic RED work used the project environment and the focused files:

```text
/home/bernt-popp/miniforge3/envs/vntyper/bin/pytest -q -m unit \
  tests/unit/test_calibration_artifacts.py \
  tests/unit/test_calibration_custody.py \
  tests/unit/test_calibration_workflow.py \
  tests/unit/test_cli_calibrate.py \
  tests/unit/test_calibration_custodian_import.py \
  tests/unit/test_calibration_study_binding.py
```

The RED phase was intentionally run in slices because the two new strict modules first
failed collection as missing. Once those interfaces were introduced, the six required
categories remained semantically red:

1. ordinary extraction still produced `locked_payload.json`, locked feature/label/
   baseline values, and its local `external-custodian` evidence manifest;
2. the old locally produced extraction-shaped package was accepted for evaluation,
   while no strict external authority binding existed;
3. `KeyboardInterrupt`/`SystemExit` did not retire interrupted locked consumption and
   exclusive writes caught only `Exception`;
4. completed `False` outcomes were not explicit and the CLI could not install their
   artifacts before exiting 1;
5. a generated profile with identical candidate parameters crossed to another study;
   and
6. atomic CLI staging survived interruption because cleanup caught only `Exception`.

During final diff review, a narrower exact-byte RED was added:

```text
/home/bernt-popp/miniforge3/envs/vntyper/bin/pytest -q -m unit \
  tests/unit/test_calibration_custodian_import.py::test_custodian_header_requires_authority_to_bind_exact_validation_bytes \
  -o log_cli=false --tb=short
```

Result before the correction: `1 failed in 0.06s` because a whitespace-altered
validation-attestation byte stream could be re-checksummed while retaining only the
semantic canonical hash. The authority now must bind the exact bytes; the affected
artifact/import suite passed `24 passed in 7.50s` afterward.

A first bare system-Python RED invocation failed collection because the system
environment lacked `pandas` and `rfc8785`. It was discarded as environment evidence and
all authoritative RED/GREEN commands used the repository's `vntyper` environment.

## GREEN verification

Final verification was run after the exact-byte correction and before the implementation
commit.

- Complete calibration plus calibration CLI family:
  `/home/bernt-popp/miniforge3/envs/vntyper/bin/pytest -q -m unit tests/unit/test_calibration_*.py tests/unit/test_cli_calibrate.py -o log_cli=false --tb=short`
  -> `316 passed in 19.91s`.
- Full branch-inclusive unit coverage gate:
  `PATH=/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH make test-unit-cov`
  -> `10150 passed, 5 skipped, 7 warnings in 346.85s`; total coverage `91.96%`
  against the `86%` floor.
- Staged patch coverage:
  `PATH=/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH make patch-coverage`
  -> `87%` over 1,648 branch-wide changed executable lines, above the `80%` gate.
- Ruff format and lint on all 15 Task 4 Python files -> `15 files already formatted` and
  `All checks passed!`.
- Mypy on the eight changed production modules -> `Success: no issues found in 8 source files`.
- `git diff --cached --check` -> clean.

The five skips and seven warnings are the repository's existing conditional skips and
dependency/FutureWarning messages, not Task 4 calibration failures.

## Files

Production:

- `vntyper/scripts/calibration_artifacts.py`
- `vntyper/scripts/calibration_custody.py`
- `vntyper/scripts/calibration_workflow.py`
- `vntyper/scripts/cli_calibrate.py`
- new `vntyper/scripts/calibration_custodian_import.py`
- new `vntyper/scripts/calibration_locked_artifacts.py`
- new `vntyper/scripts/calibration_result_artifacts.py`
- new `vntyper/scripts/calibration_study_binding.py`

Tests:

- `tests/unit/test_calibration_artifacts.py`
- `tests/unit/test_calibration_custody.py`
- `tests/unit/test_calibration_workflow.py`
- `tests/unit/test_cli_calibrate.py`
- new `tests/unit/test_calibration_custodian_import.py`
- new `tests/unit/test_calibration_result_artifacts.py`
- new `tests/unit/test_calibration_study_binding.py`

## Self-review

- Confirmed the locked payload is not read by header/inventory/checksum validation; the
  precommit is durably installed before `Path.read_bytes()` can expose it.
- Confirmed ordinary extraction's locked directory contains exactly the declaration,
  commitments, and checksum files and has no former evidence assertion.
- Confirmed the authority binds exact validation bytes as well as exact locked payload
  bytes, and every decoded study/profile/run binding is checked after the one-use open.
- Confirmed a failed or interrupted exact pair cannot retry, while another profile is
  not spuriously retired.
- Confirmed failure artifacts are installed atomically before CLI status 1 and usage
  errors remain argparse status 2.
- Confirmed every new `BaseException` catch is at an atomic cleanup or custody/
  retirement boundary and re-raises.
- Confirmed no cryptographic or independent-validation claim was introduced.

## Remaining risks and limitation

The named custodian authority is a strict externally supplied assertion, not a locally
verifiable proof of organizational independence. A future trust-root/signature policy
would be a separate reviewed contract. This task intentionally neither implements nor
claims that policy. The missing independent external cohort remains outside Task 4 and
is not resolved by these local safeguards.
