# Troubleshooting

Search this page for the first line of your error. A failed run exits 1 and its last line
is `VNtyper failed: <cause>. ... Help: <this page>`. Known causes also print a `Fix:` line.

## Where VNtyper looks for its files

**Every tool and reference path in the shipped `vntyper/config.json` is relative to the
directory you run `vntyper` from.** The defaults expect:

```text
<working directory>/
  vntyper/dependencies/kestrel/kestrel.jar, kanalyze.jar   # ship in the source checkout
  reference/                                               # vntyper install-references -d reference
```

So, for a conda or source installation:

1. run `vntyper install-references -d reference` once, **from the repository root**;
2. run every `vntyper pipeline` from that same directory (`cd /path/to/VNtyper` at the top
   of a job script);
3. after cloning a new version into a new directory, install the references there too.

The Docker image already has this layout.

## "Kestrel cannot run because required files are missing"

VNtyper checks this before reading any input, so the run stops in about a second. The
message lists each missing file, the absolute path it tried, and the working directory.
Cause: one of the three rules above was broken. Install into `reference/` in the checkout
you run from, or `cd` into the checkout that has it.

## "Kestrel reported N error(s) ... but exited 0"

The pinned Kestrel 1.0.1 exits 0 even after a fatal error, so VNtyper reads its log
(`<output>/kestrel/kestrel_kmer_<k>.log`) and quotes the error:

| Kestrel error | Fix |
| --- | --- |
| `Error reading reference sequence(s)` | Motif FASTA is empty or corrupt: `vntyper install-references -d reference` |
| `Cannot open indexed k-mer count (IKC) file`, `Error setting k-mer counts` | Check free disk space and write permission in the output directory |
| `File not found while variant writer` | The output directory was removed or is not writable |

These are never reported as a negative. Some of them leave an empty VCF behind, and an
empty VCF would otherwise read as "no variant found"
([#338](https://github.com/hassansaei/VNtyper/issues/338)).

## "Kestrel produced no usable VCF for any configured k-mer size"

Kestrel wrote no readable VCF for any k-mer size and logged no error. Read the logs the
message lists. **On versions before 2.0.38**, this error on *every* sample almost always
means the motif reference is missing. Confirm with
`grep ERROR <output>/kestrel/kestrel_kmer_20.log`, then see
[Where VNtyper looks for its files](#where-vntyper-looks-for-its-files).

## "unrecognized arguments: --reference-fasta"

`vntyper pipeline --reference-fasta` exists from **v2.0.10**. Upgrade, and always spell
options in full: older versions silently expanded `--reference` to `--reference-assembly`.

## CRAM: "Critical command failed: samtools view ..."

Decoding a CRAM needs the exact FASTA it was written against, with the same contig names and
sequences. A chr1-only FASTA or a different GRCh38 build does not work:

```bash
vntyper pipeline --cram sample.cram \
    --reference-fasta /ref/Homo_sapiens_assembly38.fasta \
    --reference-assembly hg38 -o results/sample/
```

If the FASTA lacks contigs named in the CRAM header, `pipeline.log` warns
`does not cover header contigs`. See [Input Formats](input-formats.md).

## Does a sample that fails coverage QC still get a report?

Yes. `Coverage QC: FAIL` is a finding: genotyping continues and the report shows the failed
metric. A hard error is different. A step did not run, so no result is written, rather than
showing a result that was never computed.

## The report has no IGV alignment view

**This is expected for every sample in which Kestrel called no variant.** It does not mean
that a step failed or that the installation is incomplete. The alignment view is drawn
around the position of a Kestrel call: `kestrel/output.bed` names that position, and it
is written only when a variant passes the final filter. A sample with no call has no
position to show, so it has no view. This has been the behaviour since the view was
introduced. Samples with a Kestrel call get the view.

The report's alignment panel names the reason. In `pipeline.log` it is one INFO line:

```text
INFO - No alignment view for this sample: Kestrel called no variant, and the view is drawn
around the position of a Kestrel call. This is expected for a sample without a call; ...
```

Versions up to 2.0.40 logged the same normal state as
`WARNING - BED file does not exist or not provided. Skipping IGV report generation.`
That warning was misleading and needs no action.

Three states are not normal, and each is logged as a WARNING that names what to do:

| Log starts with | Meaning | Fix |
| --- | --- | --- |
| `No alignment view although Kestrel called a variant` | The call's `kestrel/output.bed` is missing | Re-run the sample, or `vntyper report --bed-file <file>` |
| `No alignment view: Kestrel called no variant in this run, but ... exists` | An earlier run into the same `--output-dir` left its `output.bed`; it was not used | Use a fresh `--output-dir` per sample |
| `No alignment view: this run's summary holds no readable Kestrel result` | The Kestrel step was not recorded | Check the Kestrel section of the report and `pipeline.log` |

`--report-igv off` also produces no view. The report says it was switched off.

## Log messages that are normal on a successful run

A successful run ends with `Pipeline finished successfully.` and exits 0. These INFO lines
describe a sample in which Kestrel called no variant, or a stage the run did not need:

| Message | Why |
| --- | --- |
| `Kestrel called no variant: no candidate passed the final filter.` | No candidate passed the Kestrel filters, so the empty-result placeholder is written |
| `Kestrel called no variant: its VCF holds no insertion or deletion.` | The same, one step earlier |
| `Kestrel called no variant for this sample; the report states a negative Kestrel result.` | The report reading that placeholder |
| `No alignment view for this sample: Kestrel called no variant ...` | See [above](#the-report-has-no-igv-alignment-view) |
| `No fastp quality metrics in this report: fastp runs only on FASTQ input ...` | BAM and CRAM input are not run through fastp |
| `adVNTR module not included. Skipping adVNTR genotyping.` | adVNTR runs only with `--extra-modules advntr` |

`WARNING - Output directory ... already holds N entries from an earlier run (...)` names
what the `--output-dir` already contained. This run's files overwrite those with the same
name. It is harmless when you meant to re-run a sample, and one fresh directory per sample
avoids it.

Up to 2.0.40, a sample with no Kestrel call logged five of these states as WARNINGs:
`Final processed DataFrame is empty`, `No insertion/deletion variants found`,
`No Kestrel data found in pipeline summary`, `BED file does not exist or not provided` and
`fastp output file not found`. Every run also logged `Output directory ... is non-empty;
prior results may be overwritten.`, even into a fresh directory, because the run's own
`pipeline.log` was counted. None of these needs action.

## Running many samples on a cluster

Run one sample per job, so one failure cannot stop the rest. Try a single sample
interactively first, because setup errors are identical for every sample:

```bash
#!/bin/bash
#SBATCH --array=1-3500%50
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

source activate vntyper
cd /path/to/VNtyper            # config.json paths are relative to this directory

CRAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" crams.txt)
SAMPLE=$(basename "$CRAM" .cram)
vntyper pipeline --cram "$CRAM" \
    --reference-fasta /ref/Homo_sapiens_assembly38.fasta --reference-assembly hg38 \
    --output-dir "/results/vntyper/$SAMPLE" --threads 4 --fast-mode
echo -e "$SAMPLE\t$?" >> /results/vntyper/exit_status.tsv
```

Then aggregate the sample directories with [`vntyper cohort`](cohort-analysis.md).

## Keeping references outside the checkout

`--config-path` replaces the whole config, so start from a copy. When a config is given,
`install-references` writes absolute reference paths into it:

```bash
cd /path/to/VNtyper
cp vntyper/config.json /shared/vntyper/config.json
vntyper --config-path /shared/vntyper/config.json install-references -d /shared/vntyper/reference
vntyper --config-path /shared/vntyper/config.json pipeline ...
```

`tools.kestrel` and `tools.kanalyze` stay relative. Keep running from the repository root,
or set those two keys to absolute paths as well.
