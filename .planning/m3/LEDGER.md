# SDD ledger — plan: .planning/m3/PLAN.md

Branch: fix/milestone-3-web-cohort-integrity
Base: b46da80 (main)

Deviation from superpowers:subagent-driven-development, recorded deliberately:
the skill says never dispatch implementers in parallel. The user's instruction is
explicit ("Dispatch A/B/C/D as 4 concurrent subagents in one message"), and the
conflict mechanism the skill guards against — shared .git/index — is removed by
REVIEW.md finding 15: subagents run no git commands. File sets were proved
disjoint twice (my own audit + the codex review). The integrator commits.

## Progress
Task D1: complete (e331cc7, review: codex finding 1 applied, 79 tests pass)
Task A1+A3: complete (3f5c9b9) — both templates, tripwire test. Fingerprint verified unmoved.
Task A2: complete (e670f60) — js_json_literal, ensure_ascii explicit per review finding 12.

ENVIRONMENT TRAP FOUND DURING INTEGRATION: the default shell runs miniforge
*base* (plotly 6.5.2), not the `vntyper` env (plotly 6.9.0, the conda pin).
test_cohort_summary_oracle.py::test_the_cohort_report_matches_its_recorded_fingerprint
FAILS in base and PASSES in the env. Every subagent ran its tests in base.
All verification must be re-run under `conda run -n vntyper`.
