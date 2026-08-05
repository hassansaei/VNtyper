# CLAUDE.md

Project instructions for this repository live in **[AGENTS.md](AGENTS.md)**. Read it
before making any change — it is the single source of truth for setup, commands, code
style, testing, git conventions, and the repo-specific traps.

This file exists only so Claude Code picks the instructions up automatically. Keep it a
pointer: add new guidance to `AGENTS.md`, not here.

Quick reference:

- Verify with `make check-all` before proposing a PR.
- Fast feedback loop: `make format && make test-unit`, run from the repo root.
- The pipeline needs the `vntyper` conda env; `pip install -e .[dev]` alone is not enough.
