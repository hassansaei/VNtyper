# Troubleshooting

Each section below starts with the message you see, then explains the cause and the fix.
Search this page for the first line of your error.

## Where VNtyper looks for its files

Most setup errors come from this one rule. **Every tool and reference path in the shipped
`vntyper/config.json` is relative to the directory you run `vntyper` from**, not to the
config file or the installed package. The defaults expect:

```text
<working directory>/
  vntyper/dependencies/kestrel/kestrel.jar     # ships in the source checkout
  vntyper/dependencies/kestrel/kanalyze.jar    # ships in the source checkout
  reference/                                   # created by `vntyper install-references`
    All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa
    MUC1_motifs_Rev_com.fa
    alignment/chr1.hg38.fa ...
```

So a conda or source installation should:

1. run `vntyper install-references -d reference` **from the repository root**, once, and
2. run every `vntyper pipeline` from that same repository root: `cd /path/to/VNtyper` at the
   top of a job script.

Do not install the references somewhere else and leave the config unchanged. The pipeline
still looks in `./reference/`. If you need them elsewhere, see
[Keeping references outside the checkout](#keeping-references-outside-the-checkout).

The Docker image already has this layout, so none of this applies there.

## "Kestrel cannot run because required files are missing"

```text
CRITICAL - Kestrel cannot run because required files are missing. Relative paths in config.json
are resolved against the current working directory (/home/user/jobs):
  - reference_data.muc1_reference_vntr = 'reference/All_Pairwise_...fa' -> /home/user/jobs/reference/All_Pairwise_...fa (not found)
```

VNtyper checks this before it reads any input, so the run stops in under a second. The
message lists each missing file, the absolute path it tried, and the working directory it
resolved against. Usual causes:

- `install-references` was never run in this checkout, or ran with `-d` pointing somewhere
  other than `reference/`.
- The job script runs `vntyper` from a directory other than the repository root.
- A new version was cloned into a new directory, and the references stayed with the old one.

Fix: install the references into `reference/` in the checkout you run from, or `cd` into the
checkout that has them.

## "Kestrel reported N error(s) for k-mer size K but exited 0"

```text
ERROR - Kestrel reported 1 error(s) for k-mer size 20 but exited 0, so this attempt is treated
as failed rather than as a result:
  23:33:46 [KestrelRunner] ERROR e.g.kestrel.runner.KestrelRunner - Error reading reference sequence(s): ...
Full log: results/sample/kestrel/kestrel_kmer_20.log
```

The pinned Kestrel (1.0.1) exits with status 0 even after a fatal error, so VNtyper reads its
log instead. The quoted line is Kestrel's own error. The most common one,
`Error reading reference sequence(s)`, means the motif FASTA exists but is empty, truncated,
or unreadable. Reinstall it with `vntyper install-references -d reference`. Other ERROR
lines, such as `Cannot open indexed k-mer count (IKC) file`, point to a disk-space or
temporary-directory problem in the output directory.

This is never reported as a negative result. A Kestrel run that failed to load its
reference or its k-mer counts would otherwise produce an empty VCF, and an empty VCF looks
like "no variant found" ([#338](https://github.com/hassansaei/VNtyper/issues/338)).

## "Kestrel produced no usable VCF for any configured k-mer size"

Kestrel exited 0 without writing a VCF, or wrote one that has no valid header, for every
configured k-mer size, and logged no error. The message lists each attempt's log under
`<output>/kestrel/`, so read those first.

Versions **before 2.0.38** did not read the Kestrel log. There, this message on *every*
sample almost always means the motif reference is missing (see the two sections above).
Check with:

```bash
grep ERROR <output>/kestrel/kestrel_kmer_20.log
```

The run stops instead of writing a result on purpose. Writing no result, and so showing
the sample as negative, would hide a genotyper that never ran
([#212](https://github.com/hassansaei/VNtyper/issues/212)).

## "unrecognized arguments: --reference-fasta"

```text
vntyper: error: unrecognized arguments: --reference-fasta /ref/Homo_sapiens_assembly38.fasta
```

`vntyper pipeline --reference-fasta` was added in **v2.0.10**. Earlier versions have the
option only on `vntyper report`, where it means something else. Upgrade. Always spell
options in full: on those older versions argparse silently expanded `--reference` to
`--reference-assembly`, the only option it could match.

## CRAM: "Critical command failed: samtools view ..."

A CRAM stores reads as differences from the reference it was written against, so decoding
needs **that exact FASTA**, with the same contig names and sequences. Pass it explicitly:

```bash
vntyper pipeline --cram sample.cram \
    --reference-fasta /ref/Homo_sapiens_assembly38.fasta \
    --reference-assembly hg38 -o results/sample/
```

A chr1-only FASTA or a different GRCh38 build (for example UCSC `hg38.fa` for a CRAM written
against the GATK/Broad `Homo_sapiens_assembly38.fasta`) does not work. When the FASTA lacks
contigs named in the CRAM header, VNtyper logs a warning naming them. Look for
`does not cover header contigs` in `pipeline.log`. A local
reference path recorded in the CRAM header is used only when it lies inside the CRAM's own
directory. See [Input Formats](input-formats.md).

## Does a sample that fails coverage QC still get a report?

Yes. `Coverage QC: FAIL` is a finding about the sample. Genotyping continues, and the report
shows the failed metric so the call can be read in that light. A **hard error** is
different: it means a step did not run, and VNtyper refuses to present a step that did not
run as a result. The sample gets no `summary_report.html`, `pipeline.log` explains why, and
the process exits 1.

## Running many samples on a cluster

One failed sample must not stop the rest. Run one sample per job (a SLURM array), or keep
the loop going and record the exit status:

```bash
#!/bin/bash
#SBATCH --array=1-3500%50
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

source activate vntyper
cd /path/to/VNtyper            # config.json paths are relative to this directory

CRAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" crams.txt)
SAMPLE=$(basename "$CRAM" .cram)

vntyper pipeline \
    --cram "$CRAM" \
    --reference-fasta /ref/Homo_sapiens_assembly38.fasta \
    --reference-assembly hg38 \
    --output-dir "/results/vntyper/$SAMPLE" \
    --threads 4 --fast-mode
echo -e "$SAMPLE\t$?" >> /results/vntyper/exit_status.tsv
```

Then aggregate the sample directories with [`vntyper cohort`](cohort-analysis.md). Before
submitting the array, run a single sample interactively. The setup errors above are the
same for every sample and show up in about a second.

## Keeping references outside the checkout

`install-references` writes absolute paths into a config when you pass one with the global
`--config-path`. `--config-path` replaces the whole config, so start from a copy of the
shipped file:

```bash
cd /path/to/VNtyper
cp vntyper/config.json /shared/vntyper/config.json
vntyper --config-path /shared/vntyper/config.json install-references -d /shared/vntyper/reference
vntyper --config-path /shared/vntyper/config.json pipeline ...
```

The Kestrel JAR paths (`tools.kestrel`, `tools.kanalyze`) stay relative. Either keep running
from the repository root, or edit those two keys to absolute paths as well.
