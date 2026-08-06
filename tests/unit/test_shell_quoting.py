"""
Pin the shell quoting of caller-controlled paths at the sites `run_command` still builds by hand.

`command_builders.py` quotes every interpolated path with :func:`~vntyper.scripts.command_builders.quote_path`,
but several call sites were never migrated to it. `run_command` uses ``shell=True`` deliberately (it needs
process substitution for the CRAM unmapped-read path), so an unquoted operand containing a space is split
by Bash into two, and one containing a metacharacter is executed.

The quoted sites - samtools quickcheck in `utils.validate_bam_file`, and the SAM->BAM conversion, its
index and the bcftools sort in `kestrel_genotyping` - each get a path with a space *and* a single quote
in it, then assert on ``shlex.split()`` of the command handed to `run_command`: the operand list Bash
will actually see. Asserting on the split rather than on the quoted string pins the property that
matters (one path, one operand) without pinning which escaping style ``shlex.quote`` happens to choose.
Quoting a path that did not need quoting is a no-op - ``shlex.quote`` leaves anything matching
``[\\w@%+=:,./-]*`` alone - so the ordinary-path case pins that the everyday command is byte-identical
to before.

The fifth site, `install_references.execute_aligner_index`, is **still unquoted**, and the tests below
characterise it that way on purpose. `vntyper/scripts/install_references.py` is an input to the Docker
base-image content hash, so editing it forces a base-image rebuild - which a pull request cannot
publish, turning the `Build base image` job from skipped into failed. The quoting is deferred to the
next base-image rebuild, where it lands alongside issue #193 (declaring `bcrypt` directly in
`docker/requirements-web.txt`), which was held back for exactly the same reason. Characterising rather
than deleting keeps the site covered and makes the deferral visible: when the rebuild happens, these
tests fail loudly and are inverted back into the quoted form the others already use.

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


def test_aligner_index_command_does_not_yet_quote_the_reference_path(tmp_path) -> None:
    """
    Characterise the *unquoted* reference path: quoting here is deferred to the next base-image rebuild.

    The aligner templates are `.format()`ed and run with ``shell=True``, so this path is a
    caller-controlled operand exactly like the four sites above - but
    `vntyper/scripts/install_references.py` is an input to the Docker base-image content hash, and
    editing it forces a rebuild that a pull request cannot publish. The fix therefore rides along with
    the next base-image rebuild, together with issue #193.

    Until then this is what Bash actually receives: the path goes in bare, so a directory with a space
    becomes two operands and the stray single quote leaves the command unterminated - ``shlex.split``
    refuses it, which is precisely the defect the deferred quoting closes.
    """
    ref = Path(_make(tmp_path, "ref.fa"))
    aligner_info = {"index_command": "bwa index {ref_path}", "supports_threading": False}

    with patch.object(install_references.subprocess, "run") as run:
        install_references.execute_aligner_index(ref, "bwa", aligner_info)

    command = run.call_args[0][0]
    assert command == f"bwa index {ref}", "the reference path is interpolated bare, without shlex.quote"
    with pytest.raises(ValueError):
        shlex.split(command)


def test_aligner_index_derived_index_locations_are_not_yet_quoted(tmp_path) -> None:
    """
    Characterise the *unquoted* derived locations; quoting is deferred with the rest of this site.

    `index_dir` and `index_base` are derived from the reference path, so they inherit its spaces - and
    today they too are interpolated bare. See the module docstring for why the fix waits for the next
    base-image rebuild.

    With two paths on the line the two stray single quotes balance each other, so unlike the single-path
    case Bash raises nothing: it accepts the command and silently mangles it, swallowing
    ``--output-directory`` inside a quoted run. Neither path survives as an operand.
    """
    ref = Path(_make(tmp_path, "ref.fa"))
    aligner_info = {
        "index_command": "dragen-os --build-hash-table true --ht-reference {ref_path} --output-directory {index_dir}",
        "index_dir_required": True,
    }

    with patch.object(install_references.subprocess, "run") as run:
        install_references.execute_aligner_index(ref, "dragmap", aligner_info)

    index_dir = ref.parent / f"{ref.stem}_dragmap_index"
    command = run.call_args[0][0]
    assert command == f"dragen-os --build-hash-table true --ht-reference {ref} --output-directory {index_dir}"

    operands = shlex.split(command)
    assert str(ref) not in operands, "the reference path does not survive as a single operand"
    assert str(index_dir) not in operands, "the index directory does not survive as a single operand"
    assert "--output-directory" not in operands, "the flag is swallowed into a quoted run"


def test_the_thread_count_is_not_quoted(tmp_path) -> None:
    """
    Only paths are ever candidates for quoting; a numeric operand must stay a bare number.

    This one uses a shell-safe reference path, so it holds both today and after the deferred quoting
    lands - ``shlex.quote`` is a no-op on such a path, and it must never touch the thread count.
    """
    ref = tmp_path / "ref.fa"
    ref.write_text("x", encoding="utf-8")
    aligner_info = {"index_command": "bwa-mem2 index -t {threads} {ref_path}", "supports_threading": True}

    with patch.object(install_references.subprocess, "run") as run:
        install_references.execute_aligner_index(ref, "bwa-mem2", aligner_info, threads=8)

    command = run.call_args[0][0]
    assert " -t 8 " in command, "the thread count must stay a bare number"
    assert shlex.split(command) == ["bwa-mem2", "index", "-t", "8", str(ref)]
