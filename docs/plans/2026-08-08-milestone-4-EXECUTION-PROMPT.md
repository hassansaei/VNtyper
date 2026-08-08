# Execution prompt — milestone 4

Copy everything between the fences into Codex, from the repo root, on branch
`fix/milestone-4-cram-input-robustness`.

---

```
Implement VNtyper milestone 4 ("CRAM and input robustness") from its committed spec and
plan. Do not redesign: both documents survived two adversarial review rounds and every
counter-intuitive decision in them is backed by a measurement recorded in spec §3.

READ FIRST, IN THIS ORDER:
1. AGENTS.md — the source of truth for setup, style, testing and the repo's traps.
2. docs/plans/2026-08-08-milestone-4-cram-input-robustness-spec.md — WHAT and WHY.
   §3 is measured evidence. Read it before writing any CRAM code.
3. docs/plans/2026-08-08-milestone-4-cram-input-robustness-plan.md — the task list.

SUPERPOWERS — use these skills, do not skip:
- superpowers:using-git-worktrees   before Wave 2 (one worktree per parallel task)
- superpowers:subagent-driven-development  to execute the plan task by task
- superpowers:test-driven-development      inside EVERY task: RED before GREEN, no exceptions
- superpowers:verification-before-completion  before any claim that something works
- superpowers:requesting-code-review       before the PR
- superpowers:finishing-a-development-branch  to land it
Also read superpowers/skills/using-superpowers/references/codex-tools.md for the harness
adaptation.

STATE: Wave 0 (#213) is DONE and merged on this branch (commits 2d81519, 312a2d0).
Start at Wave 1, Task 1. Remaining: #225, #209, #178, #165, #161.

ORDER (the plan's waves; the dependencies are real):
  Wave 1, sequential: Tasks 1-6. Nothing in Wave 2 can start until `make check-all` is
    green at the tip of Wave 1 — every later task consumes AlignmentPlan.
  Wave 2, parallel worktrees: {Task 7} · {Task 8} · {Task 9} · {Task 10}.
    Tasks 7 and 8 in the plan are ONE worktree — #209 and #178 both edit the CRAM branch
    of process_bam_to_fastq and both edit command_builders.py. #165 and #161 are genuinely
    independent.
  Wave 3: Tasks 11-12 — golden cohort, full gate, PR, release.

FIVE MEASURED FACTS THAT LOOK WRONG AND ARE NOT. If your instinct fights one of these,
re-read spec §3 rather than "fixing" it:
1. An index built into the output directory is INVISIBLE to `samtools view` (§3.9). This
   is why the preflight builds a run-local alignment view — a symlink to the input plus a
   co-located index — instead of just passing a path around. `samtools index` writes beside
   the SYMLINK, not its target, so the read-only-input invariant holds (§3.12).
2. `samtools view -P` and `-c` CANNOT be combined (§3.10). The reference probe therefore
   runs the real slice shape into `-o /dev/null` and uses the exit status as the proof.
3. `samtools view -P` DOES fetch cross-contig mates (§3.7), so the shipped chr1-only FASTA
   is not sufficient in general, and a probe that omits `-P` authorises a reference the
   real slice then rejects.
4. A resolvable `UR:` in the CRAM header silently rescues a reference candidate that should
   fail (§3.10). Any test asserting a candidate fails MUST rename the `UR:` target first.
   Two measurements were wrong until this was noticed.
5. `samtools view -f 12 <cram> '*'` loses PLACED flag-12 reads — measured, 50 of 130 on a
   purpose-built fixture — and `samtools idxstats` column 4 reports exactly those 50, for
   BAM and CRAM alike, without needing a reference (§3.13). That is why scan selection is
   `auto` and why parsing idxstats must FAIL CLOSED: an unparsable table selects `stream`,
   never `indexed`.

EXIT BAR (spec §1): a run either succeeds, or fails BEFORE any stage does work with a
message naming exactly what is missing. No input is ever silently discarded. There is no
config value whose setting permits a read to be dropped — if you find yourself adding an
escape hatch, that is the bug, not the fix. Everything is config-driven: no threshold,
contig name, reference path or scan strategy inline in Python.

REPO TRAPS THAT WILL SILENTLY COST YOU A DAY:
- Activate the `vntyper` conda env and put its bin FIRST on PATH. The base env and
  ~/.local/bin shadow it and `make check-all` will report on the wrong tools:
    source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate vntyper
    export PATH="$CONDA_PREFIX/bin:$PATH"
- Run pytest FROM THE REPO ROOT. tests/parametrization.py opens a relative path at
  collection time.
- Every new test file needs `pytestmark = pytest.mark.unit` after the imports, or CI never
  runs it.
- A test importing app.tasks or app.main must live under tests/unit/web/ — that package's
  conftest.py is what puts docker/ on sys.path.
- Never lower the coverage floor. `make patch-coverage` must be >= 80% on changed lines.
- fastq_bam_processing.py is 649 lines and kestrel_genotyping.py is 936, both at/over the
  ~650 guideline. Changes to them must be net-negative or extract (AGENTS.md rule 3).
- `--config-path` REPLACES the whole config. Read every new key with .get and its shipped
  default.
- Inner loop: `make format && make test-unit`. Before the PR: `make check-all`. If you
  touch .github/workflows/, `make check-all` is NOT enough — run `make ci-local`.

DELIVERY:
- ONE PR against main, not a stack, with Closes #213, #225, #209, #178, #165, #161.
- Patch bump only: 2.0.9 -> 2.0.10, in ALL THREE of vntyper/version.py, CITATION.cff and
  docs/about/changelog.md (they drift and a test catches it).
- NEVER push a v*.*.* tag — it publishes to PyPI immediately and irreversibly. Release with
  `gh release create v2.0.10 --target main`.
- `gh issue view` 500s on this repo. Use `gh api repos/hassansaei/VNtyper/issues/N --jq .body`.

ADVERSARIAL GATE, still owed: round 3 was started against the spec+plan and never
completed, so both carry a round-2 verdict. Run it on the FINAL DIFF before the PR and
record the result in spec §10:
  codex exec -m gpt-5.6-sol -c model_reasoning_effort="xhigh" "Red-team this diff against
  the spec and plan in docs/plans/. Find: fixes that mask rather than remove a failure
  mode, silent-data-loss paths left open, untestable acceptance criteria, and any place the
  implementation and the spec disagree. Rank HIGH/MED/LOW with file:line evidence. Be
  maximally skeptical; do not praise."
Apply superpowers:receiving-code-review — verify each claim against the code before
accepting it, rebut what is wrong WITH EVIDENCE, fix what is right. Loop until no HIGH
survives.

REPORT AT THE END, per issue: root cause, the failure mode now impossible, the test that
proves it, and the final gate verdict. Never claim tests pass without pasting the command
output.
```
