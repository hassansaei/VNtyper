"""Every artefact path the pipeline derives from a basename, enumerated and frozen.

Contract C5. ``--output-name`` is parsed, defaulted from ``config.json`` and then
**dropped**: ``run_pipeline`` has no such parameter, and ``pipeline.py``,
``generate_report.py`` and the adVNTR stage each hardcode their own basename instead.

Threading it through looks like a one-line fix and is not, because a stale or
inconsistent artefact path does not fail hard anywhere in this pipeline:

* ``summary.record_step`` records a missing file with ``md5sum=None`` and
  ``parsed_result={"error": ...}`` rather than raising;
* ``generate_report.build_kestrel_frames`` turns that into an empty frame;
* an empty Kestrel frame is the report's configured **negative** default.

A basename that moves for some artefacts and not others is therefore a silent
false-negative genotype - exactly the failure class this work exists to prevent.

So this module enumerates the whole surface first and freezes it. Three sets, and the
split between them is the finding:

1. :data:`PARAMETERISED_ON_THE_BAM_PATH` / :data:`PARAMETERISED_ON_THE_FASTQ_PATH` -
   basenames ``run_pipeline`` hands to a stage. These *could* move, because the stage
   takes an ``output_name`` argument.
2. :data:`RECONSTRUCTED` - paths ``pipeline.py`` rebuilds from a literal instead of
   consuming what the stage returned. Contract C4 names the adVNTR pair specifically.
3. :data:`HARDCODED_IN_CONSUMERS` - basenames baked into modules ``run_pipeline`` cannot
   reach at all, so no parameter threaded through ``run_pipeline`` can move them.

Set 3 is why ``--output-name`` is refused rather than threaded; see
:mod:`vntyper.scripts.artifact_names`.
"""

import ast
import inspect
import json
from pathlib import Path
from unittest import mock

import pytest

from tests.support.pipeline_harness import artifact_paths_from_summary, run_pipeline_under_harness
from vntyper.scripts import cli_report, generate_report, pipeline
from vntyper.scripts.alignment_contract import AlignmentPlan
from vntyper.scripts.artifact_names import (
    COVERAGE_BASENAME,
    PIPELINE_BASENAME,
    pipeline_artifact_paths,
    select_best_vcf_file,
)

pytestmark = pytest.mark.unit


# --------------------------------------------------------------------------------------
# 1. Basenames handed to a stage. Every one of these stages accepts an ``output_name``.
# --------------------------------------------------------------------------------------
PARAMETERISED_ON_THE_BAM_PATH: dict[str, str] = {
    "process_bam_to_fastq": PIPELINE_BASENAME,
    "run_advntr": PIPELINE_BASENAME,
    "process_advntr_output": PIPELINE_BASENAME,
    # A *different* literal, and one nothing derives from --output-name even in
    # principle: the coverage TSV name is part of frozen contract C1.
    "calculate_vntr_coverage": COVERAGE_BASENAME,
}

PARAMETERISED_ON_THE_FASTQ_PATH: dict[str, str] = {
    "process_fastq": PIPELINE_BASENAME,
    "align_and_sort_fastq": PIPELINE_BASENAME,
    "process_bam_to_fastq": PIPELINE_BASENAME,
    "calculate_vntr_coverage": COVERAGE_BASENAME,
}

# --------------------------------------------------------------------------------------
# 2. Paths pipeline.py rebuilds itself, relative to the output directory.
# --------------------------------------------------------------------------------------
RECONSTRUCTED: set[str] = {
    "fastq_bam_processing/output_R1.fastq.gz",
    "fastq_bam_processing/output_R2.fastq.gz",
    "fastq_bam_processing/output_other.fastq.gz",
    "fastq_bam_processing/output_single.fastq.gz",
    "fastq_bam_processing/output.json",
    "fastq_bam_processing/output_sliced.bam",
    "alignment_processing/output_sorted.bam",
    "kestrel/output.vcf",
    "kestrel/output.bam",
    "kestrel/output.bed",
    "kestrel/output_indel.vcf",
    "kestrel/output_indel.vcf.gz",
    "advntr/output_adVNTR.tsv",
    "advntr/output_adVNTR.vcf",
    "advntr/output_adVNTR_result.tsv",
}

#: Artefacts whose names carry no basename at all. They never move, whatever
#: ``--output-name`` were to do, so a partial threading desynchronises them from
#: everything in :data:`RECONSTRUCTED`.
BASENAME_INDEPENDENT: set[str] = {
    "pipeline_summary.json",
    "fastq_bam_processing/pipeline_info.json",
    "coverage/coverage_summary.tsv",
    "kestrel/kestrel_result.tsv",
    "advntr/advntr_model.db",
    "advntr/cross_match_results.tsv",
    "provenance/advntr_artifact_evidence.json",
    "provenance/decision_profile.json",
}

# --------------------------------------------------------------------------------------
# 3. Basenames baked into modules run_pipeline cannot reach.
# --------------------------------------------------------------------------------------
HARDCODED_IN_CONSUMERS: dict[str, set[str]] = {
    # generate_report reads fastp's JSON by literal name; a moved basename makes
    # every fastp quality metric silently absent from the report.
    "generate_report": {"fastq_bam_processing/output.json"},
    # `vntyper report --input-dir` rediscovers the IGV inputs by literal name, and
    # the report subcommand has no --output-name flag at all, so a custom-named run
    # could never be re-reported.
    "cli_report": {"kestrel/output.bam", "kestrel/output.bed"},
}


def _relative_named_paths(harness) -> set[str]:
    """Collect every path under the output directory that a stage was handed.

    Walks the recorded call arguments rather than the filesystem: the point is what
    the pipeline *named*, not what the (stubbed) tools happened to create.

    Args:
        harness: The harness returned by ``run_pipeline_under_harness``.

    Returns:
        set[str]: POSIX-style paths relative to the output directory.
    """
    root = Path(harness.output_dir).resolve()
    found: set[str] = set()
    for stub in harness.stages.values():
        if not stub.called:
            continue
        for call in stub.call_args_list:
            for value in (*call.args, *call.kwargs.values()):
                values = value if isinstance(value, tuple) else (value,)
                for item in values:
                    candidates = (
                        (Path(item.view_path), Path(item.index_path))
                        if isinstance(item, AlignmentPlan)
                        else (Path(item),)
                        if isinstance(item, str | Path)
                        else ()
                    )
                    for candidate in candidates:
                        if not candidate.is_absolute():
                            continue
                        try:
                            found.add(candidate.resolve().relative_to(root).as_posix())
                        except ValueError:
                            continue
    return found


def _string_constants_in_source(source: str) -> set[str]:
    """Return every string literal in ``source``, **excluding docstrings**.

    ``ast.Constant`` covers docstrings too -- they are expression statements holding
    a string -- so an unfiltered walk cannot tell a filename the code *uses* from one
    a docstring merely *mentions*. ``generate_report``'s ``load_fastp_output``
    docstring says "e.g., output.json" while the production read is a single literal
    at ``generate_report.py:435``; without this filter, deleting that read left the
    test green on the strength of the prose.

    Args:
        source: Python source text.

    Returns:
        set[str]: The string constants the code evaluates.
    """
    tree = ast.parse(source)
    docstrings = {
        id(node.body[0].value)
        for node in ast.walk(tree)
        if isinstance(node, ast.Module | ast.ClassDef | ast.FunctionDef | ast.AsyncFunctionDef)
        and node.body
        and isinstance(node.body[0], ast.Expr)
        and isinstance(node.body[0].value, ast.Constant)
        and isinstance(node.body[0].value.value, str)
    }
    return {
        node.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant) and isinstance(node.value, str) and id(node) not in docstrings
    }


def _string_constants(module) -> set[str]:
    """:func:`_string_constants_in_source` applied to an imported module.

    Args:
        module: The imported module to parse.

    Returns:
        set[str]: The string constants the module's code evaluates.
    """
    return _string_constants_in_source(inspect.getsource(module))


PROSE_PROBE = '''"""Module doc naming only_in_module_doc.tsv."""


class C:
    """Class doc naming only_in_class_doc.tsv."""


def f():
    """Function doc naming only_in_func_doc.tsv."""
    return "a_real_literal.tsv"
'''


def test_the_string_constant_scan_ignores_docstrings() -> None:
    """Guard the guard: a scan that reads prose proves nothing about the code."""
    constants = _string_constants_in_source(PROSE_PROBE)

    assert constants == {"a_real_literal.tsv"}, f"a docstring leaked into the scan: {sorted(constants)}"


def test_the_scan_still_sees_a_docstrings_worth_of_text_elsewhere() -> None:
    """The filter must drop docstrings, not every long string."""
    constants = _string_constants_in_source('x = """not a docstring, an assignment"""\n')

    assert constants == {"not a docstring, an assignment"}


def test_generate_report_mentions_the_fastp_basename_in_prose_as_well_as_in_code() -> None:
    """The concrete case F4 named: without the filter this module passes on its docstring.

    ``load_fastp_output``'s docstring says "e.g., output.json" and the only production
    read spells the whole relative path. If the filter regressed, the bare basename
    would reappear in the scan and deleting line 435 would no longer fail anything.
    """
    source = inspect.getsource(generate_report)

    assert "output.json" in source, "the fastp basename left generate_report entirely"
    assert "output.json" not in _string_constants(generate_report), (
        "the bare basename is only in a docstring; if the scan sees it, the docstring filter is broken"
    )
    assert "fastq_bam_processing/output.json" in _string_constants(generate_report), (
        "the production read at generate_report.py:435 is gone; the report has lost its fastp metrics"
    )


# --------------------------------------------------------------------------------------
# The enumeration itself
# --------------------------------------------------------------------------------------


def _assert_basenames(harness, expected: dict[str, str]) -> None:
    """Assert each stage was handed the basename the contract says it gets.

    Args:
        harness: The harness returned by ``run_pipeline_under_harness``.
        expected: Stage name -> the basename it must receive.
    """
    for stage, basename in expected.items():
        stub = harness.stages[stage]
        assert stub.called, f"the pipeline never reached {stage}(); this assertion would be vacuous"
        call = stub.call_args
        actual = call.kwargs.get("output_name")
        if actual is None:
            # process_fastq, align_and_sort_fastq and the adVNTR pair take it positionally.
            positional = [value for value in call.args if isinstance(value, str)]
            assert basename in positional, f"{stage}() was not given the basename {basename!r}: {positional}"
        else:
            assert actual == basename, f"{stage}() was given basename {actual!r}, expected {basename!r}"


def test_the_bam_path_hands_every_stage_the_declared_basename(tmp_path: Path) -> None:
    """``artifact_names`` must describe the real pipeline, not an aspiration.

    Args:
        tmp_path: Pytest temporary directory.
    """
    harness = run_pipeline_under_harness(tmp_path / "out", extra_modules=["advntr"])
    _assert_basenames(harness, PARAMETERISED_ON_THE_BAM_PATH)


def test_advntr_run_snapshots_and_threads_one_verified_evidence_value(tmp_path: Path) -> None:
    from vntyper.modules.advntr.artifact_evidence import load_packaged_artifact_evidence
    from vntyper.scripts.decision_profile import load_packaged_decision_profile

    output = tmp_path / "out"
    harness = run_pipeline_under_harness(output, extra_modules=["advntr"])
    evidence = load_packaged_artifact_evidence()
    profile = load_packaged_decision_profile()
    snapshot = output / "provenance" / "advntr_artifact_evidence.json"
    profile_snapshot = output / "provenance" / "decision_profile.json"
    recorded = json.loads((output / "pipeline_summary.json").read_text(encoding="utf-8"))

    assert snapshot.read_bytes() == evidence.canonical_bytes
    assert profile_snapshot.read_bytes() == profile.canonical_bytes
    assert recorded["schema_version"] == 3
    assert recorded["advntr_evidence_digest"] == evidence.digest
    assert recorded["decision_profile_sha256"] == profile.digest
    assert recorded["decision_profile_snapshot"] == "provenance/decision_profile.json"
    parsing_evidence = harness.kwargs("process_advntr_output")["artifact_evidence"]
    reconciliation_evidence = harness.kwargs("reconcile_caller_outputs")["artifact_evidence"]
    assert parsing_evidence is reconciliation_evidence
    assert parsing_evidence.digest == evidence.digest
    from vntyper.scripts.run_configuration import resolve_run_configuration

    run_configuration = resolve_run_configuration()
    command_call = harness.kwargs("run_advntr")
    parsing_call = harness.kwargs("process_advntr_output")
    assert command_call["resolved_component"] == run_configuration.advntr
    assert command_call["runtime_component"] == run_configuration.advntr_runtime
    assert command_call["custom_context_active"] is False
    assert parsing_call["resolved_component"] == run_configuration.advntr
    assert parsing_call["nomenclature_component"] == run_configuration.nomenclature
    assert parsing_call["custom_context_active"] is False
    reconciliation_call = harness.kwargs("reconcile_caller_outputs")
    assert reconciliation_call["resolved_component"] == run_configuration.nomenclature
    cross_match_call = harness.kwargs("cross_match_variants")
    assert cross_match_call["resolved_component"] == run_configuration.cross_match
    harness.stages["load_advntr_config"].assert_not_called()


def test_pipeline_snapshots_the_supplied_explicit_context_without_reloading_package(tmp_path: Path) -> None:
    """A custom CLI context must remain the one profile named by run provenance."""
    import copy

    from vntyper.scripts.canonical_json import canonical_json_bytes, load_strict_json_object
    from vntyper.scripts.run_configuration import resolve_run_configuration

    packaged_path = Path(__file__).parents[2] / "vntyper" / "profiles" / "decision_profile.json"
    document = copy.deepcopy(load_strict_json_object(packaged_path.read_bytes()))
    document["profile_id"] = "pipeline-explicit-test"
    document["profile_revision"] = "test-1"
    document["profile_kind"] = "explicit-custom"
    inventory = document["inventory"]
    assert isinstance(inventory, dict)
    inventory["/components/kestrel/duplicate_flagging/flag_name"]["value"] = "Pipeline_Test_Duplicate"
    profile_path = tmp_path / "explicit.json"
    profile_path.write_bytes(canonical_json_bytes(document))
    run_configuration = resolve_run_configuration(profile_path)

    output = tmp_path / "out"
    harness = run_pipeline_under_harness(output, run_configuration=run_configuration, extra_modules=["advntr"])
    recorded = json.loads((output / "pipeline_summary.json").read_text(encoding="utf-8"))

    assert (output / "provenance" / "decision_profile.json").read_bytes() == (
        run_configuration.decision_profile.canonical_bytes
    )
    assert recorded["decision_profile_id"] == "pipeline-explicit-test"
    assert recorded["decision_profile_source"] == "explicit-cli"
    assert recorded["decision_profile_sha256"] == run_configuration.decision_profile.digest
    kestrel_call = harness.kwargs("run_kestrel")
    assert kestrel_call["resolved_component"] is run_configuration.kestrel
    assert kestrel_call["nomenclature_component"] is run_configuration.nomenclature
    assert kestrel_call["runtime_component"] is run_configuration.kestrel_runtime
    assert kestrel_call["custom_context_active"] is True
    parsing_call = harness.kwargs("process_advntr_output")
    assert parsing_call["resolved_component"] is run_configuration.advntr
    assert parsing_call["nomenclature_component"] is run_configuration.nomenclature
    reconciliation_call = harness.kwargs("reconcile_caller_outputs")
    assert reconciliation_call["resolved_component"] is run_configuration.nomenclature
    cross_match_call = harness.kwargs("cross_match_variants")
    assert cross_match_call["resolved_component"] is run_configuration.cross_match


def test_the_fastq_path_hands_every_stage_the_declared_basename(tmp_path: Path) -> None:
    """The FASTQ branch reaches two stages the BAM branch never does.

    Args:
        tmp_path: Pytest temporary directory.
    """
    harness = run_pipeline_under_harness(
        tmp_path / "out",
        bam=None,
        fastq1="/in/r1.fastq.gz",
        fastq2="/in/r2.fastq.gz",
    )
    _assert_basenames(harness, PARAMETERISED_ON_THE_FASTQ_PATH)


#: Every path a full BAM + adVNTR run names, frozen. This is the byte-identical
#: baseline contract C5 requires: a change here means an artefact moved, and if it
#: moved for a producer but not its consumer the result is an empty report section
#: rather than an error.
BAM_RUN_ARTEFACTS: set[str] = {
    ".",  # generate_summary_report's own output_dir argument
    "predefined_regions_hg19.bed",
    "fastq_bam_processing",
    "fastq_bam_processing/input.bam",
    "fastq_bam_processing/input.bam.bai",
    "fastq_bam_processing/output_R1.fastq.gz",
    "fastq_bam_processing/output_R2.fastq.gz",
    "fastq_bam_processing/output_other.fastq.gz",
    "fastq_bam_processing/output_single.fastq.gz",
    "fastq_bam_processing/output_sliced.bam",
    "fastq_bam_processing/pipeline_info.json",
    "coverage",
    "coverage/coverage_summary.tsv",
    "kestrel",
    "kestrel/output.vcf",
    "kestrel/output.bam",
    "kestrel/output.bed",
    "kestrel/kestrel_result.tsv",
    "advntr",
    "advntr/advntr_model.db",
    "advntr/output_adVNTR.vcf",
    "advntr/output_adVNTR_result.tsv",
    "advntr/cross_match_results.tsv",
}

#: The same, for a FASTQ run without adVNTR.
FASTQ_RUN_ARTEFACTS: set[str] = {
    ".",
    "predefined_regions_hg19.bed",
    "fastq_bam_processing",
    "fastq_bam_processing/post_alignment.bam",
    "fastq_bam_processing/post_alignment.bam.bai",
    "fastq_bam_processing/output.json",
    "fastq_bam_processing/output_R1.fastq.gz",
    "fastq_bam_processing/output_R2.fastq.gz",
    "fastq_bam_processing/output_other.fastq.gz",
    "fastq_bam_processing/output_single.fastq.gz",
    "alignment_processing",
    "alignment_processing/output_sorted.bam",
    "coverage",
    "coverage/coverage_summary.tsv",
    "kestrel",
    "kestrel/output.vcf",
    "kestrel/output.bam",
    "kestrel/output.bed",
    "kestrel/kestrel_result.tsv",
}


def test_the_bam_run_artefact_set_is_unchanged(tmp_path: Path) -> None:
    """The frozen enumeration for BAM input plus adVNTR, compared all at once.

    Args:
        tmp_path: Pytest temporary directory.
    """
    harness = run_pipeline_under_harness(tmp_path / "out", extra_modules=["advntr"])
    named = _relative_named_paths(harness) | artifact_paths_from_summary(harness.output_dir)
    assert named == BAM_RUN_ARTEFACTS, (
        f"artefact paths drifted; new: {sorted(named - BAM_RUN_ARTEFACTS)}; gone: {sorted(BAM_RUN_ARTEFACTS - named)}"
    )


def test_the_fastq_run_artefact_set_is_unchanged(tmp_path: Path) -> None:
    """The frozen enumeration for FASTQ input, which names two more artefacts.

    Args:
        tmp_path: Pytest temporary directory.
    """
    out = tmp_path / "out"
    harness = run_pipeline_under_harness(
        out,
        bam=None,
        fastq1="/in/r1.fastq.gz",
        fastq2="/in/r2.fastq.gz",
    )
    named = _relative_named_paths(harness) | artifact_paths_from_summary(out)
    assert named == FASTQ_RUN_ARTEFACTS, (
        f"artefact paths drifted; new: {sorted(named - FASTQ_RUN_ARTEFACTS)}; "
        f"gone: {sorted(FASTQ_RUN_ARTEFACTS - named)}"
    )


def test_the_declaration_is_exactly_the_enumerated_surface(tmp_path: Path) -> None:
    """The declaration partitions into the two documented halves and nothing else.

    :data:`RECONSTRUCTED` is what moves with the basename; :data:`BASENAME_INDEPENDENT`
    is what never does. Anything outside both is an artefact this analysis missed.

    Args:
        tmp_path: Pytest temporary directory.
    """
    declared = {Path(path).relative_to(tmp_path).as_posix() for path in pipeline_artifact_paths(tmp_path).values()}
    # Produced inside process_bam_to_fastq and fastp; the pipeline never names them,
    # so they are absent from the run-based sets above but are still artefacts.
    produced_but_never_named = {
        "fastq_bam_processing/output_other.fastq.gz",
        "fastq_bam_processing/output_single.fastq.gz",
        "fastq_bam_processing/output.html",
    }
    assert declared == RECONSTRUCTED | BASENAME_INDEPENDENT | produced_but_never_named, (
        f"unclassified: {sorted(declared - RECONSTRUCTED - BASENAME_INDEPENDENT - produced_but_never_named)}; "
        f"classified but not declared: "
        f"{sorted((RECONSTRUCTED | BASENAME_INDEPENDENT | produced_but_never_named) - declared)}"
    )


def test_every_path_a_run_names_is_declared(tmp_path: Path) -> None:
    """``pipeline_artifact_paths`` must not miss anything a run actually produces.

    This is what binds :mod:`vntyper.scripts.artifact_names` to ``pipeline.py``: an
    artefact named by the pipeline but absent from the declaration is one that a
    future basename change would silently leave behind.

    Args:
        tmp_path: Pytest temporary directory.
    """
    declared = set(pipeline_artifact_paths(tmp_path / "out").values())
    declared_relative = {Path(path).relative_to(tmp_path / "out").as_posix() for path in declared}

    named: set[str] = set()
    for kwargs in (
        {"extra_modules": ["advntr"]},
        {"bam": None, "fastq1": "/in/r1.fastq.gz", "fastq2": "/in/r2.fastq.gz"},
    ):
        run_dir = tmp_path / f"out{len(named)}"
        harness = run_pipeline_under_harness(run_dir, **kwargs)
        named |= {
            path
            for path in _relative_named_paths(harness) | artifact_paths_from_summary(run_dir)
            if PIPELINE_BASENAME in path or COVERAGE_BASENAME in path
        }

    assert named, "no basename-derived paths were named; this test would be vacuous"
    # The stage directories are handed to stages as output directories, not as
    # artefacts, and "coverage" is a directory whose name happens to be a basename.
    stage_dirs = {"coverage", "kestrel", "advntr", "fastq_bam_processing", "alignment_processing"}
    undeclared = named - declared_relative - stage_dirs - {"in.bam"}
    assert not undeclared, f"named by a run but absent from pipeline_artifact_paths: {sorted(undeclared)}"


def test_every_declared_path_moves_together_under_a_custom_basename(tmp_path: Path) -> None:
    """The declaration's own consistency: change the basename, everything moves.

    This is the property a threaded ``--output-name`` would need, and the reason it
    cannot be threaded: the declaration can satisfy it, but three consumers named in
    :data:`HARDCODED_IN_CONSUMERS` cannot.

    Args:
        tmp_path: Pytest temporary directory.
    """
    default = pipeline_artifact_paths(tmp_path, basename=PIPELINE_BASENAME)
    moved = pipeline_artifact_paths(tmp_path, basename="mysample")

    assert set(default) == set(moved), "the declaration lost or gained a key when the basename changed"
    unmoved = {key for key, path in default.items() if moved[key] == path}
    # Only the basename-independent artefacts may stay put.
    assert unmoved == {
        "coverage_summary",
        "kestrel_result",
        "advntr_model",
        "cross_match",
        "pipeline_info",
        "pipeline_summary",
        "advntr_evidence",
        "decision_profile",
    }, f"these did not move with the basename: {sorted(unmoved)}"


def test_the_kestrel_igv_inputs_are_rebuilt_from_a_literal(tmp_path: Path) -> None:
    """``pipeline.py`` names Kestrel's outputs itself rather than consuming them.

    ``run_kestrel`` returns nothing; ``pipeline.py`` reconstructs ``output.bam``,
    ``output.bed`` and ``output_indel.vcf`` from the same literal Kestrel hardcodes
    internally. Moving one without the other silently drops the IGV panel.

    Args:
        tmp_path: Pytest temporary directory.
    """
    out = tmp_path / "out"
    (out / "kestrel").mkdir(parents=True)
    (out / "kestrel" / "output_indel.vcf").touch()

    harness = run_pipeline_under_harness(out)
    kwargs = harness.kwargs("generate_summary_report")

    assert kwargs["bam_file"] == str(out / "kestrel" / f"{PIPELINE_BASENAME}.bam")
    assert kwargs["bed_file"] == str(out / "kestrel" / f"{PIPELINE_BASENAME}.bed")
    assert kwargs["vcf_file"] == str(out / "kestrel" / f"{PIPELINE_BASENAME}_indel.vcf")


def test_the_compressed_vcf_wins_when_both_exist(tmp_path: Path) -> None:
    """The other half of ``_select_best_vcf_file``'s literal pair.

    Args:
        tmp_path: Pytest temporary directory.
    """
    kestrel_dir = tmp_path / "kestrel"
    kestrel_dir.mkdir()
    (kestrel_dir / f"{PIPELINE_BASENAME}_indel.vcf").touch()
    (kestrel_dir / f"{PIPELINE_BASENAME}_indel.vcf.gz").touch()

    assert select_best_vcf_file(str(kestrel_dir)) == str(kestrel_dir / f"{PIPELINE_BASENAME}_indel.vcf.gz")


@pytest.mark.parametrize("input_kind", ["bam", "cram", "fastq"])
def test_archive_protects_original_cli_inputs_not_the_routed_tuple(input_kind: str, tmp_path: Path) -> None:
    """Post-conversion reassignment must not replace operator-owned archive protections."""
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    original = {name: inputs / name for name in ("sample.bam", "sample.cram", "R1.fastq.gz", "R2.fastq.gz")}
    reference = inputs / "reference.fa"
    bed = inputs / "regions.bed"
    bwa = inputs / "bwa.fa"
    for path in (*original.values(), reference, bed, bwa):
        path.touch()
    bam: str | None = None
    cram: str | None = None
    fastq1: str | None = None
    fastq2: str | None = None
    if input_kind == "bam":
        bam = str(original["sample.bam"])
    elif input_kind == "cram":
        cram = str(original["sample.cram"])
    else:
        fastq1 = str(original["R1.fastq.gz"])
        fastq2 = str(original["R2.fastq.gz"])

    with mock.patch.object(pipeline, "create_safe_archive", return_value=str(tmp_path / "out.zip")) as archive:
        run_pipeline_under_harness(
            tmp_path / "out",
            archive_results=True,
            reference_fasta=reference,
            bed_file=bed,
            bwa_reference=str(bwa),
            bam=bam,
            cram=cram,
            fastq1=fastq1,
            fastq2=fastq2,
        )

    protected = archive.call_args.kwargs["protected_paths"]
    expected_inputs = {
        "bam": {str(original["sample.bam"])},
        "cram": {str(original["sample.cram"])},
        "fastq": {str(original["R1.fastq.gz"]), str(original["R2.fastq.gz"])},
    }[input_kind]
    assert set(protected) == expected_inputs | {reference, bed, str(bwa)}
    assert all(not isinstance(path, tuple) for path in protected)


# --------------------------------------------------------------------------------------
# Why the parameter is refused rather than threaded
# --------------------------------------------------------------------------------------


def test_the_configured_default_basename_is_the_one_the_pipeline_uses() -> None:
    """``config.json`` must not advertise a basename the pipeline ignores.

    ``default_values.output_name`` was ``"processed"`` while every artefact was named
    ``output``. Because the CLI dropped the value before calling ``run_pipeline``, the
    disagreement was invisible - and it is exactly what makes threading the parameter
    unsafe: threading the configured default would have moved every path on an
    ordinary run.
    """
    with open(Path(pipeline.__file__).parent.parent / "config.json", encoding="utf-8") as handle:
        config = json.load(handle)
    assert config["default_values"]["output_name"] == PIPELINE_BASENAME


@pytest.mark.parametrize("module_name", sorted(HARDCODED_IN_CONSUMERS))
def test_the_consumers_hardcode_the_basename_out_of_run_pipelines_reach(module_name: str) -> None:
    """Pin the coupling that makes a threaded ``--output-name`` unverifiable.

    These modules are not called with a basename and take no parameter for one. Any
    ``--output-name`` threaded through ``run_pipeline`` would move the producer and
    leave these readers behind, which is a silently empty report section.

    Args:
        module_name: ``generate_report`` or ``cli_report``.
    """
    module = {"generate_report": generate_report, "cli_report": cli_report}[module_name]
    constants = _string_constants(module)
    assert constants, f"parsed no string constants out of {module_name}; this test would be vacuous"

    hits = {path for path in HARDCODED_IN_CONSUMERS[module_name] if any(part in constants for part in _parts(path))}
    assert hits == HARDCODED_IN_CONSUMERS[module_name], (
        f"{module_name} no longer hardcodes {sorted(HARDCODED_IN_CONSUMERS[module_name] - hits)}; "
        "if it now takes the basename as a parameter, --output-name may finally be threadable"
    )


def _parts(path: str) -> set[str]:
    """Return the spellings a module might use for ``path``.

    Args:
        path: A POSIX-style relative artefact path.

    Returns:
        set[str]: The whole path and its final component.
    """
    return {path, Path(path).name}
