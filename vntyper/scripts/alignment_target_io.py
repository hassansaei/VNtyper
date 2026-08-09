"""Safe installation and ownership validation for generated pipeline files."""

from __future__ import annotations

import logging
import os
import stat
import tempfile
from contextlib import suppress
from pathlib import Path

from vntyper.scripts.alignment_contract import AlignmentPlan, index_candidate_names
from vntyper.scripts.reference_resolution import configured_reference_candidates

logger = logging.getLogger(__name__)

BWA_INDEX_EXTENSIONS = (".amb", ".ann", ".bwt", ".pac", ".sa")


def _reject(message: str) -> None:
    logger.error(message)
    raise ValueError(message)


def _absolute(path: str | Path) -> Path:
    return Path(os.path.abspath(path))


def _same_file(left: Path, right: Path) -> bool:
    try:
        return os.path.samefile(left, right)
    except OSError:
        return False


def _protected_alignment_paths(input_path: str | Path, file_format: str) -> tuple[Path, ...]:
    paths = [_absolute(input_path)]
    paths.extend(
        _absolute(candidate)
        for candidate in index_candidate_names(str(input_path), file_format)
        if os.path.lexists(candidate)
    )
    return tuple(paths)


def _protected_path_variants(paths: tuple[str | Path, ...]) -> tuple[Path, ...]:
    variants: list[Path] = []
    for path in paths:
        absolute = _absolute(path)
        variants.extend((absolute, absolute.resolve(strict=False)))
    return tuple(dict.fromkeys(variants))


def bwa_index_paths(reference: str | Path, config: dict) -> tuple[Path, ...]:
    """Return the BWA sidecars used by alignment and pipeline safety.

    Args:
        reference: BWA reference FASTA.
        config: Pipeline configuration containing optional index suffixes.

    Returns:
        The canonical sidecars plus any configured extras adjacent to ``reference``.
    """
    configured = config.get("tool_params", {}).get("bwa_index_extensions", ())
    extensions = tuple(dict.fromkeys((*BWA_INDEX_EXTENSIONS, *configured)))
    reference_path = Path(reference)
    return tuple(reference_path.with_name(reference_path.name + extension) for extension in extensions)


def reference_index_paths(reference: str | Path) -> tuple[Path, ...]:
    """Return the FASTA sidecars protected with a CRAM reference.

    Args:
        reference: Reference FASTA path.

    Returns:
        The samtools/htslib FASTA index, compressed FASTA index, and sequence
        dictionary paths derived from ``reference``.
    """
    reference_path = Path(reference)
    return (
        Path(f"{reference_path}.fai"),
        Path(f"{reference_path}.gzi"),
        reference_path.with_suffix(".dict"),
    )


def alignment_operator_paths(
    input_path: str | Path,
    file_format: str,
    bed_file: str | Path | None,
    reference_fasta: str | Path | None,
    config: dict,
    reference_assembly: str,
) -> tuple[str | Path, ...]:
    """Enumerate alignment-mode paths that remain operator owned.

    Args:
        input_path: Selected BAM or CRAM.
        file_format: Alignment format, ``bam`` or ``cram``.
        bed_file: Optional operator-provided BED.
        reference_fasta: Optional explicit CRAM reference.
        config: Pipeline configuration containing CRAM reference candidates.
        reference_assembly: Assembly label selecting exact or family references.

    Returns:
        The alignment, existing index candidates, BED, references, and reference
        sidecars protected by the pre-write boundary.
    """
    protected: list[str | Path] = list(_protected_alignment_paths(input_path, file_format))
    if bed_file is not None:
        protected.append(bed_file)
    references: list[str | Path] = [reference_fasta] if reference_fasta is not None else []
    if file_format == "cram":
        references.extend(
            value
            for _source, value in configured_reference_candidates(config, reference_assembly)
            if isinstance(value, str)
        )
    for reference_input in references:
        reference = Path(reference_input)
        protected.extend((reference, *reference_index_paths(reference)))
    return tuple(dict.fromkeys(protected))


def fastq_operator_paths(
    fastq_1: str | Path,
    fastq_2: str | Path | None,
    bed_file: str | Path | None,
    bwa_reference: str | Path,
    config: dict,
) -> tuple[str | Path, ...]:
    """Enumerate direct-FASTQ paths that remain operator owned.

    Args:
        fastq_1: First FASTQ input.
        fastq_2: Optional second FASTQ input.
        bed_file: Optional operator-provided BED.
        bwa_reference: BWA reference FASTA.
        config: Pipeline configuration containing optional BWA index suffixes.

    Returns:
        FASTQs, BED, BWA reference, and every canonical or configured sidecar.
    """
    reference = Path(bwa_reference)
    protected: list[str | Path] = [fastq_1, reference, *bwa_index_paths(reference, config)]
    if fastq_2 is not None:
        protected.append(fastq_2)
    if bed_file is not None:
        protected.append(bed_file)
    return tuple(dict.fromkeys(protected))


def _alignment_plan_protected_paths(
    plan: AlignmentPlan,
    additional_paths: tuple[str | Path, ...],
) -> tuple[Path, ...]:
    paths: list[Path] = []
    for alignment_path in (plan.input_path, plan.view_path):
        alignment = _absolute(alignment_path)
        paths.extend((alignment, alignment.resolve(strict=False)))
        for source in (alignment, alignment.resolve(strict=False)):
            paths.extend(_absolute(candidate) for candidate in index_candidate_names(str(source), plan.file_format))
    index = _absolute(plan.index_path)
    paths.extend((index, index.resolve(strict=False)))
    operator_paths = additional_paths
    if plan.reference_path is not None:
        operator_paths = (*operator_paths, plan.reference_path)
    paths.extend(_protected_path_variants(operator_paths))
    return tuple(dict.fromkeys(paths))


def _validate_owned_destination(
    destination: Path,
    protected_paths: tuple[Path, ...],
    *,
    description: str = "derived alignment-conversion",
    allow_symlink: bool = False,
) -> None:
    destination_absolute = _absolute(destination)
    if os.path.lexists(destination):
        destination_stat = os.lstat(destination)
        if stat.S_ISLNK(destination_stat.st_mode):
            if allow_symlink:
                return
            _reject(f"Unsafe {description} destination is a symlink: {destination}")
        if not stat.S_ISREG(destination_stat.st_mode) and not stat.S_ISLNK(destination_stat.st_mode):
            _reject(f"Unsafe {description} destination is not a regular file: {destination}")
    else:
        destination_stat = None

    for protected in protected_paths:
        if destination_absolute == protected or _same_file(destination, protected):
            _reject(f"Unsafe {description} destination aliases protected input: {destination}")

    if destination_stat is not None and destination_stat.st_nlink > 1:
        _reject(f"Unsafe {description} destination has multiple hard links: {destination}")


def _alignment_conversion_destinations(output: str | Path, output_name: str) -> tuple[Path, ...]:
    output_root = Path(output)
    sliced_bam = output_root / f"{output_name}_sliced.bam"
    destinations = [
        sliced_bam,
        output_root / f"{output_name}_unmapped.bam",
        output_root / f"{output_name}_sliced_unmapped.bam",
        output_root / f"{output_name}_R1.fastq.gz",
        output_root / f"{output_name}_R2.fastq.gz",
        output_root / f"{output_name}_other.fastq.gz",
        output_root / f"{output_name}_single.fastq.gz",
        output_root / f"{output_name}_slice.log",
        output_root / f"{output_name}_filter.log",
        output_root / f"{output_name}_merge.log",
        output_root / f"{output_name}_index.log",
        output_root / f"{output_name}_sort_fastq.log",
    ]
    destinations.extend(Path(candidate) for candidate in index_candidate_names(str(sliced_bam), "bam"))
    return tuple(destinations)


def _fastq_processing_destinations(output: str | Path, output_name: str, *, paired: bool) -> tuple[Path, ...]:
    output_root = Path(output)
    destinations = [
        output_root / f"{output_name}_R1.fastq.gz",
        output_root / f"{output_name}.html",
        output_root / f"{output_name}.json",
        output_root / f"{output_name}_fastp.log",
    ]
    if paired:
        destinations.append(output_root / f"{output_name}_R2.fastq.gz")
    return tuple(destinations)


def _is_within(path: Path, root: Path) -> bool:
    return path == root or root in path.parents


def _validate_operator_inputs_outside_output(
    output_root: Path,
    operator_inputs: tuple[str | Path, ...],
) -> tuple[Path, ...]:
    protected_paths = _protected_path_variants(operator_inputs)
    root_variants = _protected_path_variants((output_root,))
    for protected in protected_paths:
        if any(_is_within(protected, root) for root in root_variants):
            _reject(f"Operator-owned input is inside pipeline output root: {protected}")
    return protected_paths


def _validate_existing_output_tree(
    output_root: Path,
    protected_paths: tuple[Path, ...],
    *,
    allowed_file_symlink_paths: set[Path],
) -> None:
    if not os.path.lexists(output_root):
        return
    if not output_root.is_dir():
        _reject(f"Pipeline output root must be a directory: {output_root}")

    pending = [output_root]
    while pending:
        directory = pending.pop()
        with os.scandir(directory) as entries:
            for entry in entries:
                path = Path(entry.path)
                if entry.is_symlink():
                    if _absolute(path) in allowed_file_symlink_paths:
                        continue
                    if entry.is_dir(follow_symlinks=True):
                        _reject(f"Unsafe pipeline output tree contains a symlink directory: {path}")
                    _reject(f"Unsafe pipeline output tree contains an unsafe file symlink: {path}")
                if entry.is_dir(follow_symlinks=False):
                    pending.append(path)
                    continue
                if not entry.is_file(follow_symlinks=False):
                    _reject(f"Unsafe pipeline output tree contains a non-regular entry: {path}")
                if any(_same_file(path, protected) for protected in protected_paths):
                    _reject(f"Pipeline output-tree entry aliases protected input: {path}")
                if entry.stat(follow_symlinks=False).st_nlink > 1:
                    _reject(f"Unsafe pipeline output tree contains a file with multiple hard links: {path}")


def protect_alignment_inputs(
    output: str | Path,
    input_path: str | Path,
    file_format: str,
    bed_file: str | Path | None,
    reference_fasta: str | Path | None,
    config: dict,
    reference_assembly: str,
) -> None:
    """Protect alignment-mode inputs before any pipeline-owned write.

    The output tree is checked before validation logs, target preparation, preflight,
    or stage-directory creation. Existing single-link regular rerun artifacts are
    accepted, as is the exact run-local preflight view symlink. Its builder performs
    the next ownership check and atomically replaces that entry regardless of its
    pre-existing target.

    Args:
        output: Pipeline output root.
        input_path: Operator-owned BAM or CRAM.
        file_format: Alignment format, ``bam`` or ``cram``.
        bed_file: Operator-provided BED, when present.
        reference_fasta: Operator-provided CRAM reference, when present.
        config: Pipeline configuration containing CRAM reference candidates.
        reference_assembly: Assembly label selecting exact or family reference keys.

    Raises:
        ValueError: If an operator input or existing output-tree entry is unsafe.
    """
    output_root = Path(output)
    validate_alignment_output_root(output_root, input_path, file_format)
    protected_inputs = alignment_operator_paths(
        input_path,
        file_format,
        bed_file,
        reference_fasta,
        config,
        reference_assembly,
    )
    protected_paths = _validate_operator_inputs_outside_output(output_root, protected_inputs)
    preflight_view = output_root / "fastq_bam_processing" / f"input.{file_format}"
    _validate_existing_output_tree(
        output_root,
        protected_paths,
        allowed_file_symlink_paths={_absolute(preflight_view)},
    )


def validate_alignment_conversion_destinations(
    output: str | Path,
    output_name: str,
    plan: AlignmentPlan,
    *,
    protected_inputs: tuple[str | Path, ...] = (),
) -> None:
    """Validate every file path owned by alignment-to-FASTQ conversion.

    Existing single-link regular files are replaceable run artifacts. Symlinks,
    non-regular entries, multiply linked files, and aliases of the proven plan are
    rejected without modifying any destination.

    Args:
        output: Conversion-stage output directory.
        output_name: Base name used for every derived artifact.
        plan: Proven alignment inputs and index that destinations must not alias.
        protected_inputs: Additional operator-owned inputs, such as a provided BED.

    Raises:
        ValueError: If an output root or derived destination is unsafe.
    """
    validate_alignment_output_root(output, plan.input_path, plan.file_format)
    protected_paths = _alignment_plan_protected_paths(plan, protected_inputs)
    for destination in _alignment_conversion_destinations(output, output_name):
        _validate_owned_destination(destination, protected_paths)


def validate_fastq_processing_destinations(
    output: str | Path,
    output_name: str,
    fastq_1: str | Path,
    fastq_2: str | Path | None,
) -> None:
    """Reject fastp destinations that could overwrite operator-owned FASTQs.

    Existing single-link regular files remain valid rerun artifacts. Direct path
    collisions and filesystem aliases are rejected before fastp is launched.

    Args:
        output: FASTQ-processing output directory.
        output_name: Base name used by fastp outputs and reports.
        fastq_1: Operator-owned first FASTQ input.
        fastq_2: Optional operator-owned second FASTQ input.

    Raises:
        ValueError: If a fastp destination is unsafe or aliases an input FASTQ.
    """
    protected_inputs: tuple[str | Path, ...] = (fastq_1,)
    if fastq_2 is not None:
        protected_inputs = (*protected_inputs, fastq_2)
    protected_paths = _protected_path_variants(protected_inputs)
    for destination in _fastq_processing_destinations(output, output_name, paired=fastq_2 is not None):
        _validate_owned_destination(destination, protected_paths, description="FASTQ quality-control")


def validate_fastq_pipeline_destinations(
    output: str | Path,
    fastq_1: str | Path,
    fastq_2: str | Path | None,
    bed_file: str | Path | None,
    bwa_reference: str | Path,
    config: dict,
) -> None:
    """Validate the full direct-FASTQ destination set against every input.

    The validation runs once before stage-directory creation or any external tool,
    covering destinations owned by fastp, BWA, alignment preflight, and post-alignment
    conversion. Existing regular files must have one link. The exact post-alignment
    view symlink remains replaceable regardless of its pre-existing target because
    its builder revalidates and atomically replaces it before stage reads begin.

    Args:
        output: Pipeline output root.
        fastq_1: Operator-owned first FASTQ.
        fastq_2: Optional operator-owned second FASTQ.
        bed_file: Operator-provided BED target, or ``None`` for a generated target.
        bwa_reference: BWA reference FASTA.
        config: Pipeline configuration containing optional BWA index suffixes.

    Raises:
        ValueError: If any destructive destination is unsafe or aliases an input.
    """
    output_root = Path(output)
    fastq_output = output_root / "fastq_bam_processing"
    alignment_output = output_root / "alignment_processing"
    sorted_bam = alignment_output / "output_sorted.bam"
    post_alignment = fastq_output / "post_alignment.bam"
    protected_inputs = fastq_operator_paths(fastq_1, fastq_2, bed_file, bwa_reference, config)
    protected_paths = _validate_operator_inputs_outside_output(output_root, protected_inputs)
    _validate_existing_output_tree(
        output_root,
        protected_paths,
        allowed_file_symlink_paths={_absolute(post_alignment)},
    )
    destinations = [
        *_fastq_processing_destinations(fastq_output, "output", paired=fastq_2 is not None),
        sorted_bam,
        sorted_bam.with_suffix(".bam.bai"),
        alignment_output / "output_alignment.log",
        alignment_output / "output_index.log",
        post_alignment,
        *(Path(candidate) for candidate in index_candidate_names(str(post_alignment), "bam")),
        fastq_output / "post_alignment_index.log",
        fastq_output / "post_alignment_idxstats.log",
        fastq_output / "post_alignment_alignment_probe.log",
        *_alignment_conversion_destinations(fastq_output, "output"),
    ]
    for destination in dict.fromkeys(destinations):
        _validate_owned_destination(
            destination,
            protected_paths,
            description="direct-FASTQ pipeline",
            allow_symlink=destination == post_alignment,
        )


def remove_validated_slice_indexes(output: str | Path, output_name: str) -> None:
    """Remove stale slice indexes only after destination validation succeeds.

    All candidate entries are rechecked before the first unlink so a concurrent
    unsafe replacement cannot cause partial cleanup. Paths are derived directly
    beneath the conversion output directory.

    Args:
        output: Conversion-stage output directory.
        output_name: Base name used for the sliced BAM.

    Raises:
        ValueError: If a candidate escaped the output or became unsafe.
    """
    output_root = _absolute(output)
    sliced_bam = output_root / f"{output_name}_sliced.bam"
    candidates = tuple(Path(candidate) for candidate in index_candidate_names(str(sliced_bam), "bam"))
    for candidate in candidates:
        if _absolute(candidate).parent != output_root:
            _reject(f"Slice index cleanup escaped the conversion output: {candidate}")
        _validate_owned_destination(candidate, ())
    for candidate in candidates:
        if os.path.lexists(candidate):
            candidate.unlink()


def validate_alignment_output_root(
    output_root: str | Path,
    input_path: str | Path,
    file_format: str,
) -> None:
    """Reject an output root that could modify an alignment's input tree.

    Args:
        output_root: Intended pipeline or generated-target output directory.
        input_path: Protected BAM/CRAM input.
        file_format: Alignment format used to enumerate protected indexes.

    Raises:
        ValueError: If the root is unsafe or aliases protected input state.
    """
    root = Path(output_root)
    root_absolute = _absolute(root)
    root_variants = (root_absolute, root_absolute.resolve(strict=False))
    logical_input_tree = _absolute(input_path).parent
    resolved_input_tree = _absolute(input_path).resolve(strict=False).parent
    for root_variant in root_variants:
        for input_tree in (logical_input_tree, resolved_input_tree):
            if root_variant == input_tree or input_tree in root_variant.parents:
                _reject(f"Alignment output root must stay outside the patient input tree: {root}")

    if os.path.lexists(root) and not root.is_dir():
        _reject(f"Alignment output root must be a directory: {root}")
    for protected in _protected_alignment_paths(input_path, file_format):
        if root_absolute == protected or _same_file(root, protected):
            _reject(f"Alignment output root aliases protected alignment input or index: {root}")


def _validate_generated_destination(
    destination: Path,
    input_path: str | Path | None,
    file_format: str | None,
) -> None:
    destination_absolute = _absolute(destination)
    if os.path.lexists(destination):
        if destination.is_symlink():
            _reject(f"Generated BED destination must not be a symlink: {destination}")
        if not stat.S_ISREG(os.lstat(destination).st_mode):
            _reject(f"Generated BED destination must be a regular file: {destination}")

    if input_path is None or file_format is None:
        return

    validate_alignment_output_root(destination.parent, input_path, file_format)
    protected_paths = _protected_alignment_paths(input_path, file_format)
    for protected in protected_paths:
        if destination_absolute == protected or _same_file(destination, protected):
            _reject(f"Generated BED destination aliases protected alignment input or index: {destination}")


def install_generated_bed(
    destination: str | Path,
    bed_text: str,
    *,
    input_path: str | Path | None,
    file_format: str | None,
) -> None:
    """Atomically install generated BED text without following a destination link.

    Args:
        destination: Final generated BED path.
        bed_text: Fully formatted BED contents.
        input_path: Protected BAM/CRAM path, or ``None`` for FASTQ input.
        file_format: Protected alignment format, or ``None`` for FASTQ input.

    Raises:
        OSError: If the temporary write or atomic replacement fails.
        ValueError: If the output could mutate the input tree or an unsafe entry.
    """
    final_path = Path(destination)
    _validate_generated_destination(final_path, input_path, file_format)
    previous_mode = stat.S_IMODE(os.lstat(final_path).st_mode) if os.path.lexists(final_path) else None
    descriptor, temporary_name = tempfile.mkstemp(
        dir=final_path.parent,
        prefix=f".{final_path.name}.",
        suffix=".tmp",
    )
    temporary_path = Path(temporary_name)
    try:
        if previous_mode is not None:
            os.fchmod(descriptor, previous_mode)
        with os.fdopen(descriptor, "w", encoding="utf-8") as temporary_file:
            descriptor = -1
            temporary_file.write(bed_text)
            temporary_file.flush()
            os.fsync(temporary_file.fileno())
        os.replace(temporary_path, final_path)
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        with suppress(OSError):
            os.unlink(temporary_path)
