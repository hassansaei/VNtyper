"""
Pin the shell quoting of caller-controlled paths at the sites `run_command` still builds by hand.

`command_builders.py` quotes every interpolated path with :func:`~vntyper.scripts.command_builders.quote_path`,
but several call sites were never migrated to it. `run_command` uses ``shell=True`` deliberately - every
command it is given is shell syntax (pipes, ``&&``, redirects) and the ``set -o pipefail`` prefix is a
non-POSIX ``set`` option that pinning ``/bin/bash`` guarantees - so an unquoted operand containing a space
is split by Bash into two, and one containing a metacharacter is executed.

(That rationale used to be stated as "it needs process substitution for the CRAM unmapped-read path". The
CRAM branch no longer uses process substitution - it is a plain pipe, see
:func:`~vntyper.scripts.command_builders.build_cram_unmapped_filter_command` - but ``shell=True`` and
``executable="/bin/bash"`` are still required, and so is the quoting this module pins.)

All five sites are quoted: samtools quickcheck in `utils.validate_bam_file`, the SAM->BAM conversion, its
index and the bcftools sort in `kestrel_genotyping`, and the aligner indexing command in
`install_references.execute_aligner_index`. Each test below gets a path with a space *and* a single
quote in it, then asserts on ``shlex.split()`` of the command handed to `run_command` (or, for the last
one, `subprocess.run`): the operand list Bash will actually see. Asserting on the split rather than on
the quoted string pins the property that matters (one path, one operand) without pinning which escaping
style ``shlex.quote`` happens to choose. Quoting a path that did not need quoting is a no-op -
``shlex.quote`` leaves anything matching ``[\\w@%+=:,./-]*`` alone - so the ordinary-path case pins that
the everyday command is byte-identical to before.

The fifth site, `install_references.execute_aligner_index`, was quoted last, on purpose:
`vntyper/scripts/install_references.py` is an input to the Docker base-image content hash, so editing it
forces a base-image rebuild - which a pull request cannot publish, turning the `Build base image` job
from skipped into failed. The quoting was deferred to the next base-image rebuild, where it landed
alongside issue #193 (declaring `bcrypt` directly in `docker/requirements-web.txt`), which was held back
for exactly the same reason. It quotes with `install_references._quote`, a three-line copy of
:func:`~vntyper.scripts.command_builders.quote_path` rather than an import of it, because
`install_references.py` may import nothing from `vntyper` at module scope (see the import comment at the
top of that file); `TestQuote` in `tests/unit/test_install_references_derivations.py` pins that the two
behave identically.

Deliberately *not* quoted in any case: `config["tools"]` entries and the aligner ``index_command``
templates. Those are command prefixes such as ``mamba run -n envadvntr advntr``, which must stay
multiple words.
"""

import os
import shlex
from pathlib import Path
from unittest.mock import patch

import pytest

from vntyper.scripts import install_references, kestrel_genotyping, utils

pytestmark = pytest.mark.unit

#: A path fragment that is split by Bash unless quoted, and whose quote would unbalance the command.
HOSTILE = "od d's dir"


def _make(tmp_path: Path, name: str) -> str:
    """Create ``<tmp_path>/<HOSTILE>/<name>`` and return it as a string."""
    directory = tmp_path / HOSTILE
    directory.mkdir(exist_ok=True)
    target = directory / name
    target.write_text("x", encoding="utf-8")
    return str(target)


def test_bam_quickcheck_quotes_the_alignment_path(tmp_path) -> None:
    """A BAM under a directory with a space must reach samtools as one operand."""
    bam = _make(tmp_path, "sample.bam")

    with patch.object(utils, "run_command", return_value=True) as run:
        utils.validate_bam_file(bam)

    command = run.call_args[0][0]
    assert shlex.split(command) == ["samtools", "quickcheck", "-v", bam]


def test_bam_quickcheck_leaves_an_ordinary_path_alone(tmp_path) -> None:
    """Quoting is a no-op for shell-safe paths, so the usual command is byte-identical."""
    bam = str(tmp_path / "sample.bam")
    Path(bam).write_text("x", encoding="utf-8")

    with patch.object(utils, "run_command", return_value=True) as run:
        utils.validate_bam_file(bam)

    assert run.call_args[0][0] == f"samtools quickcheck -v {bam}"


def test_sam_to_bam_quotes_both_operands_and_the_redirect(tmp_path) -> None:
    """The redirect target is a shell word too: an unquoted space there writes two files."""
    sam = _make(tmp_path, "output.sam")
    output_dir = str(tmp_path / HOSTILE)

    with patch.object(kestrel_genotyping, "run_command", return_value=True) as run:
        kestrel_genotyping.convert_sam_to_bam_and_index(sam, output_dir)

    view, index = run.call_args_list[0][0][0], run.call_args_list[1][0][0]
    bam = os.path.join(output_dir, "output.bam")
    assert shlex.split(view) == ["samtools", "view", "-Sb", sam, ">", bam]
    assert shlex.split(index) == ["samtools", "index", bam]


def test_sam_to_bam_uses_configured_samtools_and_threads_in_both_commands(tmp_path) -> None:
    """Both converter commands use the configured executable and non-default threads."""
    sam = str(tmp_path / "output.sam")
    Path(sam).write_text("@HD\tVN:1.6\n", encoding="utf-8")

    with patch.object(kestrel_genotyping, "run_command", return_value=True) as run:
        kestrel_genotyping.convert_sam_to_bam_and_index(
            sam,
            str(tmp_path),
            samtools_path="/opt/vntyper/bin/samtools",
            threads=7,
        )

    bam = str(tmp_path / "output.bam")
    assert run.call_args_list[0].args[0] == f"/opt/vntyper/bin/samtools view -Sb -@ 7 {sam} > {bam}"
    assert run.call_args_list[1].args[0] == f"/opt/vntyper/bin/samtools index -@ 7 {bam}"


def test_bcftools_sort_quotes_the_vcf_paths(tmp_path) -> None:
    """`-o` takes a single operand; an unquoted space makes bcftools see a stray argument."""
    vcf = _make(tmp_path, "out.vcf")
    gz = _make(tmp_path, "out.vcf.gz")

    with (
        patch.object(kestrel_genotyping.shutil, "which", return_value="/usr/bin/bcftools"),
        patch.object(kestrel_genotyping, "run_command", return_value=True) as run,
    ):
        kestrel_genotyping._try_compress_vcf_with_bcftools(vcf, gz, str(tmp_path))

    assert shlex.split(run.call_args[0][0]) == ["bcftools", "sort", vcf, "-o", gz, "-W", "-O", "z"]


def test_aligner_index_command_quotes_the_reference_path(tmp_path) -> None:
    """A reference path with a space and a stray quote must reach the aligner as one operand.

    The aligner templates are `.format()`ed and run with ``shell=True``, so this path is a
    caller-controlled operand exactly like the four sites above, quoted with
    `install_references._quote` before the template is formatted.
    """
    ref = Path(_make(tmp_path, "ref.fa"))
    aligner_info = {"index_command": "bwa index {ref_path}", "supports_threading": False}

    with patch.object(install_references.subprocess, "run") as run:
        install_references.execute_aligner_index(ref, "bwa", aligner_info)

    command = run.call_args[0][0]
    assert shlex.split(command) == ["bwa", "index", str(ref)]


def test_aligner_index_derived_index_locations_are_quoted(tmp_path) -> None:
    """`index_dir` and `index_base` are derived from the reference path, so they inherit its spaces
    and are quoted the same way it is - each survives as its own operand."""
    ref = Path(_make(tmp_path, "ref.fa"))
    aligner_info = {
        "index_command": "dragen-os --build-hash-table true --ht-reference {ref_path} --output-directory {index_dir}",
        "index_dir_required": True,
    }

    with patch.object(install_references.subprocess, "run") as run:
        install_references.execute_aligner_index(ref, "dragmap", aligner_info)

    index_dir = ref.parent / f"{ref.stem}_dragmap_index"
    command = run.call_args[0][0]
    assert shlex.split(command) == [
        "dragen-os",
        "--build-hash-table",
        "true",
        "--ht-reference",
        str(ref),
        "--output-directory",
        str(index_dir),
    ]


def test_the_thread_count_is_not_quoted(tmp_path) -> None:
    """
    Only paths are ever candidates for quoting; a numeric operand must stay a bare number.

    Uses a reference path with a space and a stray quote so the same command line proves both things
    at once: the path reaches bwa-mem2 as one operand, and the thread count next to it stays bare.
    """
    ref = Path(_make(tmp_path, "ref.fa"))
    aligner_info = {"index_command": "bwa-mem2 index -t {threads} {ref_path}", "supports_threading": True}

    with patch.object(install_references.subprocess, "run") as run:
        install_references.execute_aligner_index(ref, "bwa-mem2", aligner_info, threads=8)

    command = run.call_args[0][0]
    assert " -t 8 " in command, "the thread count must stay a bare number"
    assert shlex.split(command) == ["bwa-mem2", "index", "-t", "8", str(ref)]
