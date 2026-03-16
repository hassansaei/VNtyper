# CLI Reference

Overview of VNtyper's command-line interface.

**Usage:** `vntyper [global-options] <command> [command-options]`

Global options can be placed before or after the subcommand.

## Global Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-v, --version` | flag | — | Show version and exit |
| `-l, --log-level` | choice | `INFO` | Logging level: `DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL` |
| `-f, --log-file` | path | (stdout) | Log output file path. For `pipeline` and `cohort` commands, automatically set to `<output-dir>/pipeline.log` or `<output-dir>/cohort.log` when `--output-dir` is provided |
| `--config-path` | path | (bundled config) | Path to a custom `config.json`. If not provided, the default bundled configuration is used |

## Commands

| Command | Description |
|---------|-------------|
| [`pipeline`](pipeline.md) | Run the full genotyping pipeline |
| [`install-references`](install-references.md) | Download and set up reference files |
| [`report`](report.md) | Generate HTML summary report |
| [`cohort`](cohort.md) | Aggregate multi-sample results |
| [`online`](online.md) | Submit analysis to vntyper.org |
