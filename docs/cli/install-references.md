# install-references

Download and set up the reference files required by VNtyper 2.

## Synopsis

```
vntyper [global-options] install-references
    -d <dir>
    [--skip-indexing]
    [-t <int>]
    [--aligners <aligner> ...]
    [--references <reference> ...]
```

## Arguments

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-d, --output-dir` | path | (required) | Directory where references will be installed |
| `--skip-indexing` | flag | off | Skip the indexing step |
| `-t, --threads` | int | `4` | Number of threads to use for indexing |
| `--aligners` | list | `bwa` only | Specific aligners to use (e.g., `bwa`, `bwa-mem2`, `minimap2`). If not specified, only BWA is used |
| `--references` | list | `hg19 hg38` | Specific references to process (e.g., `hg19`, `hg38`, `GRCh37`, `GRCh38`). Default: UCSC references only (hg19, hg38) |

## Examples

Install references with default settings:

```bash
vntyper install-references -d /path/to/refs
```

Install with multiple aligners and 8 threads:

```bash
vntyper install-references -d /path/to/refs --aligners bwa bwa-mem2 -t 8
```

Install specific references without indexing:

```bash
vntyper install-references -d /path/to/refs --references GRCh37 GRCh38 --skip-indexing
```
