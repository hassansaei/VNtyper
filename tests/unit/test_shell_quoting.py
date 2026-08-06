"""
Pin the shell quoting of caller-controlled paths at the sites `run_command` still builds by hand.

`command_builders.py` quotes every interpolated path with :func:`~vntyper.scripts.command_builders.quote_path`,
but four call sites were never migrated to it. `run_command` uses ``shell=True`` deliberately (it needs
process substitution for the CRAM unmapped-read path), so an unquoted operand containing a space is split
by Bash into two, and one containing a metacharacter is executed.

Every test here builds a path with a space *and* a single quote in it, then asserts on
``shlex.split()`` of the command handed to `run_command` - the operand list Bash will actually see.
Asserting on the split rather than on the quoted string pins the property that matters (one path,
one operand) without pinning which escaping style ``shlex.quote`` happens to choose. Quoting a path
that did not need quoting is a no-op - ``shlex.quote`` leaves anything matching ``[\\w@%+=:,./-]*``
alone - so the ordinary-path case pins that the everyday command is byte-identical to before.

Deliberately *not* quoted: `config["tools"]` entries and the aligner ``index_command`` templates. Those are
command prefixes such as ``mamba run -n envadvntr advntr``, which must stay multiple words.
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
    """
    The aligner templates are `.format()`ed and run with ``shell=True``.

    The *template* stays unquoted - it is a command prefix - but every path substituted into it is
    a caller-controlled operand.
    """
    ref = Path(_make(tmp_path, "ref.fa"))
    aligner_info = {"index_command": "bwa index {ref_path}", "supports_threading": False}

    with patch.object(install_references.subprocess, "run") as run:
        install_references.execute_aligner_index(ref, "bwa", aligner_info)

    assert shlex.split(run.call_args[0][0]) == ["bwa", "index", str(ref)]


def test_aligner_index_command_quotes_the_derived_index_locations(tmp_path) -> None:
    """`index_dir` and `index_base` are derived from the reference path, so they inherit its spaces."""
    ref = Path(_make(tmp_path, "ref.fa"))
    aligner_info = {
        "index_command": "dragen-os --build-hash-table true --ht-reference {ref_path} --output-directory {index_dir}",
        "index_dir_required": True,
    }

    with patch.object(install_references.subprocess, "run") as run:
        install_references.execute_aligner_index(ref, "dragmap", aligner_info)

    index_dir = ref.parent / f"{ref.stem}_dragmap_index"
    assert shlex.split(run.call_args[0][0]) == [
        "dragen-os",
        "--build-hash-table",
        "true",
        "--ht-reference",
        str(ref),
        "--output-directory",
        str(index_dir),
    ]


def test_the_thread_count_is_not_quoted(tmp_path) -> None:
    """Only paths are quoted; a numeric operand must stay a bare number."""
    ref = Path(_make(tmp_path, "ref.fa"))
    aligner_info = {"index_command": "bwa-mem2 index -t {threads} {ref_path}", "supports_threading": True}

    with patch.object(install_references.subprocess, "run") as run:
        install_references.execute_aligner_index(ref, "bwa-mem2", aligner_info, threads=8)

    command = run.call_args[0][0]
    assert " -t 8 " in command, "the thread count must stay a bare number"
    assert shlex.split(command) == ["bwa-mem2", "index", "-t", "8", str(ref)]
