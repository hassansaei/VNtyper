# #238 — the generated reference sidecar is replaced by a link to itself

CRAM preflight fails under Docker with `Too many levels of symbolic links` (ELOOP) while
opening the run-local reference index. This is one defect with a two-line fix and a
defence-in-depth layer, not an environment problem.

## 1. The one root cause

`ReferenceBinding.bind_generated_sidecars()` retains an htslib-generated `.fai`/`.gzi` by
opening it and then **replacing that same pathname** with a symlink to
`/proc/<pid>/fd/<n>` — a descriptor whose own recorded path *is* the pathname just
replaced. The entry becomes self-referential:

```
.input_reference_1/reference.fa.fai  ->  /proc/7/fd/12
/proc/7/fd/12                        ->  .../.input_reference_1/reference.fa.fai (deleted)
```

Nothing but Linux's dentry-jump implementation of procfs magic links keeps that
resolvable. A runtime that resolves `/proc/<pid>/fd/<n>` by re-walking its link text
closes the cycle and every consumer that opens the sidecar gets ELOOP.

This is the **only** binding whose source and destination are the same path.
`AlignmentBinding` and every construction-time `_InodeView` point at an operator-owned
path elsewhere on disk, so no cycle can form. That asymmetry is exactly what the
reporter's log shows.

## 2. Measured evidence

### 2.1 The reporter's log is explained line by line

| Log line | Candidate | Sidecar origin | Outcome |
| --- | --- | --- | --- |
| `refs_load_fai … reference.fa.fai: Too many levels of symbolic links` | 1 (`--reference-fasta`) | **generated** by htslib — the operator FASTA had no `.fai` | ELOOP |
| `[W::cram_get_ref] … ref 'chr2' not present` | 3 (shipped `chr1.hg38.fa`) | **shipped** `.fai`, bound at construction | reference loaded; fails only on genuine chr2 absence |

Candidate 3 proves magic links resolve fine in that runtime. Only the candidate whose
sidecar went through `bind_generated_sidecars()` failed.

The 11-second gap between `07:15:36` and the first error is htslib indexing the 3.2 GB
FASTA during the target probe — the step that creates the sidecar in the first place.

The two `input.cram: Too many levels of symbolic links` lines are stale-`errno` noise:
`hts_open_format` prints `strerror(errno)` after the reference load already failed.
Candidate 3 opens the same `input.cram` one second later and decodes slices from it.

### 2.2 The cycle, observed in a real run

`vntyper:latest` (2.0.12), reporter's command shape, hg38 CRAM, `--reference-fasta`
pointing at a FASTA with **no** `.fai` beside it — so htslib generates one. Watched from
inside the container while the pipeline ran:

```
SIDECAR ENTRY : …/fastq_bam_processing/.input_reference_1/reference.fa.fai
  -> points to: /proc/22/fd/7
  -> that FD's own path: …/.input_reference_1/reference.fa.fai (deleted)
```

The entry points at a descriptor whose own path is that entry. The retained inode has
`nlink=0`; the only name for it is a link back to itself. On Linux the ` (deleted)`
marker and dentry jumping hide the cycle. After the fix, the same watch on the same
image with only the three changed files swapped in reports:

```
SIDECAR ENTRY : …/.input_reference_1/reference.fa.fai
  -> TYPE: regular file, links=1, size=23
```

and the pipeline finishes successfully.

### 2.3 What is *not* the cause — ruled out by measurement

| Hypothesis | Test | Result |
| --- | --- | --- |
| Native Linux Docker, same double bind-mount | reporter's exact command, `vntyper:latest` (2.0.12), hg38 CRAM | **pipeline completes** |
| QEMU user-mode emulation (`--platform`) | procfd view probe under `linux/arm64` on x86_64 | resolves fine |
| FUSE-backed bind mount | passthrough FUSE + procfd view probe | resolves fine |

### 2.4 The fallback could not have saved the run

With `-v .:/opt/vntyper/input -v .:/opt/vntyper/output` the two mounts report the **same
`st_dev`** yet `os.link()` across them fails:

```
same st_dev across the two bind mounts: True
hardlink fallback: FAILED [18] Invalid cross-device link
```

So `(st_dev, st_ino)` is not a valid test for "same filesystem": the procfd view is the
only strategy available in the documented Docker layout, and when it fails there is
nothing behind it.

## 3. The fix

**3.1 Stop self-replacing the generated sidecar.** `_replace_generated_entry()` keeps the
regular entry and records its identity — the branch the code already takes when the proc
link cannot be installed. Replacing buys nothing: after the swap the name is a symlink,
so anything that could replace the name could replace the link too, and the sidecar lives
in a `0700` directory this binding creates and removes.

**3.2 Prove every installed view through its own path.** Both binding classes verified
the procfs *source* and never the installed *destination*, so a view no external tool
can open was still published. `consumer_reachable_identity()` now opens the destination
the way a child process would and requires `(st_dev, st_ino)` to equal the bound
descriptor identity. A view that fails the proof is logged with its `errno` and the
hardlink strategy is used instead; if that also fails, the existing `RuntimeError` names
the path and the reason rather than leaving htslib to emit ELOOP text later.

This layer does not fix #238 on its own — §2.4 shows the hardlink cannot be created in
the reporter's mount layout — but it converts any future "installed view is unusable"
into one actionable line at install time.

## 4. The second defect the same command exposes

`validate_alignment_output_root()` compared pathnames only — the lexical form and
`resolve()` — never `os.path.samefile`, which lives in the same module and is used
elsewhere in it. A directory mounted at two container paths has two different names and
one inode, so the patient-input-tree guard silently did not fire and the run wrote its
output tree inside the directory holding the patient alignment.

The bypass is an accident of the mount layout, not a policy: the identical layout with a
**single** mount (`-v .:/data`, `--cram /data/x.cram -o /data/out`) was already rejected,
because there the lexical comparison sees `/data` among `/data/out`'s parents.

### 4.1 Measured

```
same st_dev across the two bind mounts: True
os.path.samefile(input_mount, output_mount): True
container-path guard would compare: /opt/vntyper/output/run vs /opt/vntyper/input
resolve(output_run): /opt/vntyper/output/run
resolve(input dir):  /opt/vntyper/input
```

### 4.2 The fix

`_aliased_input_ancestor()` walks the output root and its parents and compares each
against the logical and resolved input trees by inode. A hit is rejected with a message
naming both pathnames, because the whole difficulty is that they look unrelated:

```
Alignment output root must stay outside the patient input tree:
/opt/vntyper/output/NTI_test.cram_v2.0.12 lies under /opt/vntyper/output, which is the
same directory as the patient input tree /opt/vntyper/input.
Give the run separate input and output directories.
```

Verified against real bind mounts, not simulated ones: the reporter's command shape is
now rejected with exactly that message, and the documented separate-directory layout
(`-v "$PWD":/opt/vntyper/input -v "$PWD/results":/opt/vntyper/output`) completes with
`Pipeline finished successfully`.

Both the `README.md` and `docs/user-guide/docker.md` volume-mount sections now state the
requirement, show the rejected form, and show the accepted one.
