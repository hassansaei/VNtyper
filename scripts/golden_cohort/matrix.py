"""Derive the golden-cohort case matrix from what is actually in ``tests/data``.

The gate page describes the matrix in prose: "the 7 multi-reference samples at all six
assemblies plus their original hg19 subsets, and the hg38 regression guard
``example_40cf``", then five repeats without ``--fast-mode`` and three ``--extra-modules
advntr`` runs. A hardcoded list of 58 ids reproduces that only until ``tests/data``
changes, and then reproduces it wrongly and silently. So the 50 BAM-by-assembly cases are
**derived** by walking the data directory, and the three selections that are policy rather
than data - which cases repeat without fast mode, which run adVNTR, which are the
deliberate-mismatch probes - are declared here by case id and **resolved** against the
derived set. A policy naming a case the data does not contain is an error, not a silent
drop.

What the page pins and what it does not
---------------------------------------
The adVNTR selection is named on the page (``a5c1_hg19_advntr``, ``b178_hg19_advntr``,
``dfc3_hg19_advntr``). The non-fast selection is **not**: the page says only "five cases
repeat without ``--fast-mode``". Two constraints on it are recoverable, and
:data:`NON_FAST_CASE_IDS` satisfies both:

* Run 2's assembly-guard table (20 hg19 / 9 hg38 / 8 GRCh38 / 7 each GRCh37,
  hg19_ensembl, hg38_ensembl) minus the 50 derived cases and the 3 adVNTR cases leaves
  exactly 3 hg19, 1 hg38 and 1 GRCh38 non-fast cases.
* Adjudication D1 counts 10 negative cases, which resolves to ``example_7a61`` at seven
  assemblies plus ``example_40cf`` at hg38 plus **one** non-fast rerun of each - so the
  non-fast set contains exactly one ``7a61`` case and exactly one ``40cf`` case.

The individual ids are a reconstruction, not a recovery. They are declared here so the
choice is visible and can be overridden with ``--non-fast-cases`` rather than being
buried in a lost script.

Attributes:
    ASSEMBLIES: The six assemblies the cohort is provided at, in the page's order.
    NON_FAST_CASE_IDS: Which derived cases repeat without ``--fast-mode``.
    ADVNTR_CASE_IDS: Which derived cases repeat with ``--extra-modules advntr``.
    PROBE_SPECS: The three deliberate-mismatch probes.
    DOCUMENTED_ASSEMBLY_COUNTS: Run 2's per-assembly case counts, used as a self-check.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

#: The six assemblies, in the order the gate page lists them.
ASSEMBLIES: tuple[str, ...] = ("hg19", "hg38", "GRCh37", "GRCh38", "hg19_ensembl", "hg38_ensembl")

#: ``example_<id>_<assembly>_subset.bam`` directly under the data directory.
SUBSET_BAM_RE = re.compile(r"^(?P<sample>example_[0-9A-Za-z]+)_(?P<assembly>.+)_subset\.bam$")

#: The leading ``example_<id>`` of a remapped BAM; its assembly comes from the directory.
REMAPPED_SAMPLE_RE = re.compile(r"^(?P<sample>example_[0-9A-Za-z]+)_")

#: See the module docstring: reconstructed, not recovered.
NON_FAST_CASE_IDS: tuple[str, ...] = (
    "7a61_hg19_subset",
    "a5c1_hg19_subset",
    "b178_hg19_subset",
    "40cf_hg38_subset",
    "dfc3_GRCh38_bwa",
)

#: Named on the gate page as ``a5c1_hg19_advntr`` / ``b178_hg19_advntr`` /
#: ``dfc3_hg19_advntr``, all at ``--advntr-max-coverage 300``.
ADVNTR_CASE_IDS: tuple[str, ...] = (
    "a5c1_hg19_subset",
    "b178_hg19_subset",
    "dfc3_hg19_subset",
)

#: ``(probe_id, base_case_id, declared_assembly, expectation)``. Two deliberate
#: mismatches, which the page records as exit 1 on both sides with only the failure point
#: moving, and one naming probe (an NCBI-named BAM declared ``hg38``) which exits 0.
PROBE_SPECS: tuple[tuple[str, str, str, str], ...] = (
    ("probe_mismatch_hg19_as_hg38", "dfc3_hg19_subset", "hg38", "nonzero"),
    ("probe_mismatch_hg38_as_hg19", "40cf_hg38_subset", "hg19", "nonzero"),
    ("probe_naming_ncbi_as_hg38", "dfc3_GRCh38_bwa", "hg38", "zero"),
)

#: Run 2's after-side assembly-guard verdict counts, which are also the per-assembly case
#: counts. Used only as a self-check on the derivation; a mismatch is reported loudly and
#: does not stop the run, because ``tests/data`` is allowed to grow.
DOCUMENTED_ASSEMBLY_COUNTS: dict[str, int] = {
    "hg19": 20,
    "hg38": 9,
    "GRCh38": 8,
    "GRCh37": 7,
    "hg19_ensembl": 7,
    "hg38_ensembl": 7,
}

#: The page's own totals, checked the same way.
DOCUMENTED_TOTALS: dict[str, int] = {"base": 50, "nonfast": 5, "advntr": 3, "total": 58, "probes": 3}


def _short(sample: str) -> str:
    """Strip the ``example_`` prefix a case id does not carry.

    Args:
        sample: The sample name as it appears in the BAM filename.

    Returns:
        str: The short form used in case ids.
    """
    return sample.removeprefix("example_")


def derive_base_cases(data_dir: Path) -> tuple[list[dict[str, Any]], list[str]]:
    """Walk ``tests/data`` and build one case per BAM-by-assembly it provides.

    Two shapes are recognised, and nothing else is: ``example_<id>_<assembly>_subset.bam``
    directly under ``data_dir``, and ``remapped/<aligner>/<assembly>/example_<id>*.bam``.
    A BAM that matches neither is logged and skipped rather than guessed at.

    Args:
        data_dir: The ``tests/data`` directory.

    Returns:
        tuple[list[dict], list[str]]: The derived cases sorted by id, and the human-readable
        derivation log describing what was found and what was skipped.

    Raises:
        ValueError: If ``data_dir`` does not exist.
    """
    if not data_dir.is_dir():
        msg = f"Test data directory not found: {data_dir}. Run `make download-test-data` first."
        logger.error(msg)
        raise ValueError(msg)

    log: list[str] = [f"data_dir={data_dir}"]
    cases: list[dict[str, Any]] = []

    for bam in sorted(data_dir.glob("*.bam")):
        match = SUBSET_BAM_RE.match(bam.name)
        if match is None:
            log.append(f"skipped (name does not match the subset pattern): {bam.name}")
            continue
        sample = match.group("sample")
        assembly = match.group("assembly")
        cases.append(
            {
                "case_id": f"{_short(sample)}_{assembly}_subset",
                "kind": "pipeline",
                "group": "base",
                "sample": sample,
                "assembly": assembly,
                "source": "subset",
                "bam": str(bam.resolve()),
                "fast_mode": True,
                "advntr": False,
                "expect_exit": "zero",
            }
        )
    log.append(f"subset BAMs: {sum(1 for c in cases if c['source'] == 'subset')}")

    remapped_root = data_dir / "remapped"
    if remapped_root.is_dir():
        for bam in sorted(remapped_root.glob("*/*/*.bam")):
            aligner = bam.parent.parent.name
            assembly = bam.parent.name
            match = REMAPPED_SAMPLE_RE.match(bam.name)
            if match is None:
                log.append(f"skipped (no example_<id> prefix): {bam.relative_to(data_dir)}")
                continue
            sample = match.group("sample")
            cases.append(
                {
                    "case_id": f"{_short(sample)}_{assembly}_{aligner}",
                    "kind": "pipeline",
                    "group": "base",
                    "sample": sample,
                    "assembly": assembly,
                    "source": f"remapped/{aligner}",
                    "bam": str(bam.resolve()),
                    "fast_mode": True,
                    "advntr": False,
                    "expect_exit": "zero",
                }
            )
        log.append(f"remapped BAMs: {sum(1 for c in cases if c['source'].startswith('remapped'))}")
    else:
        log.append(f"no remapped/ directory under {data_dir}; only subset BAMs were derived")

    cases.sort(key=lambda case: case["case_id"])
    log.append(f"derived base cases: {len(cases)}")
    return cases, log


def _resolve(case_ids: tuple[str, ...], by_id: dict[str, dict[str, Any]], what: str) -> list[dict[str, Any]]:
    """Look policy-named case ids up in the derived set, refusing any that is absent.

    Args:
        case_ids: The declared ids.
        by_id: The derived cases, keyed by id.
        what: What the policy selects, for the error message.

    Returns:
        list[dict]: The named base cases, in the declared order.

    Raises:
        ValueError: If any id is not in the derived set. Silently dropping one would
            shrink the matrix without changing a line of the page it is checked against.
    """
    missing = [case_id for case_id in case_ids if case_id not in by_id]
    if missing:
        msg = (
            f"The {what} policy names {len(missing)} case(s) that tests/data does not provide: "
            f"{', '.join(missing)}. Either the data changed or the policy is stale; fix one of them "
            "rather than running a smaller matrix than the gate page describes."
        )
        logger.error(msg)
        raise ValueError(msg)
    return [by_id[case_id] for case_id in case_ids]


def apply_policies(
    base_cases: list[dict[str, Any]],
    *,
    non_fast_ids: tuple[str, ...],
    advntr_ids: tuple[str, ...],
    advntr_max_coverage: int,
) -> tuple[list[dict[str, Any]], list[str]]:
    """Add the non-fast and adVNTR repeats to the derived base cases.

    Args:
        base_cases: The derived BAM-by-assembly cases.
        non_fast_ids: Which base cases repeat without ``--fast-mode``.
        advntr_ids: Which base cases repeat with ``--extra-modules advntr``.
        advntr_max_coverage: The value passed to ``--advntr-max-coverage``.

    Returns:
        tuple[list[dict], list[str]]: All pipeline cases, and the derivation log lines.
    """
    by_id = {case["case_id"]: case for case in base_cases}
    log: list[str] = []
    cases = list(base_cases)

    for base in _resolve(non_fast_ids, by_id, "non-fast"):
        case = dict(base)
        case.update(
            {
                "case_id": f"{_short(base['sample'])}_{base['assembly']}_nonfast",
                "group": "nonfast",
                "fast_mode": False,
                "repeat_of": base["case_id"],
            }
        )
        cases.append(case)
    log.append(f"non-fast repeats: {len(non_fast_ids)} -> {', '.join(non_fast_ids)}")

    for base in _resolve(advntr_ids, by_id, "adVNTR"):
        case = dict(base)
        case.update(
            {
                "case_id": f"{_short(base['sample'])}_{base['assembly']}_advntr",
                "group": "advntr",
                "advntr": True,
                "advntr_max_coverage": advntr_max_coverage,
                "repeat_of": base["case_id"],
            }
        )
        cases.append(case)
    log.append(f"adVNTR repeats: {len(advntr_ids)} at --advntr-max-coverage {advntr_max_coverage}")

    return cases, log


def build_probes(base_cases: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Build the three deliberate-mismatch probes from :data:`PROBE_SPECS`.

    Args:
        base_cases: The derived BAM-by-assembly cases.

    Returns:
        list[dict]: The probe cases. Probes are run and compared but are not part of the
        58-case matrix, exactly as the page counts them.
    """
    by_id = {case["case_id"]: case for case in base_cases}
    probes: list[dict[str, Any]] = []
    for probe_id, base_id, declared, expectation in PROBE_SPECS:
        base = _resolve((base_id,), by_id, f"probe {probe_id}")[0]
        probe = dict(base)
        probe.update(
            {
                "case_id": probe_id,
                "group": "probe",
                "assembly": declared,
                "declared_assembly": declared,
                "true_assembly": base["assembly"],
                "expect_exit": expectation,
                "repeat_of": base["case_id"],
            }
        )
        probes.append(probe)
    return probes


def build_cohort_cases(pipeline_case_ids: list[str]) -> list[dict[str, Any]]:
    """Build the cohort-mode cases over the per-sample output directories.

    The gate page records cohort mode as **not covered**, which is why these exist. The
    flags are read off ``vntyper/scripts/cli_parser.py``: ``-i/--input-dirs`` (or
    ``--input-file``, one of which is required), ``-o/--output-dir`` (required),
    ``--summary-file``, ``--summary-formats`` and ``--pseudonymize-samples``.

    Four cases, each covering something the others do not:

    * ``cohort_multi`` - every per-sample directory, all three export formats.
    * ``cohort_multi_pseudonymized`` - the same inputs with ``--pseudonymize-samples``,
      which is the only way ``pseudonymization_table.tsv`` is written at all.
    * ``cohort_single`` - one sample, the smallest cohort the CLI accepts as input.
    * ``cohort_empty`` - a directory holding no ``pipeline_summary.json``. The CLI cannot
      be given zero input directories (the group is ``required=True``), so this is the
      smallest empty case there is; ``aggregate_cohort`` logs an error and **returns**,
      writing no report and exiting 0.

    Args:
        pipeline_case_ids: The ids of the per-sample cases that will have run first.

    Returns:
        list[dict]: The cohort cases, in run order.
    """
    single = pipeline_case_ids[:1]
    return [
        {
            "case_id": "cohort_multi",
            "kind": "cohort",
            "group": "cohort",
            "inputs": list(pipeline_case_ids),
            "summary_formats": "csv,tsv,json",
            "pseudonymize": None,
            "expect_exit": "zero",
            "allow_missing_inputs": False,
        },
        {
            "case_id": "cohort_multi_pseudonymized",
            "kind": "cohort",
            "group": "cohort",
            "inputs": list(pipeline_case_ids),
            "summary_formats": "csv,tsv,json",
            "pseudonymize": "sample_",
            "expect_exit": "zero",
            "allow_missing_inputs": False,
        },
        {
            "case_id": "cohort_single",
            "kind": "cohort",
            "group": "cohort",
            "inputs": single,
            "summary_formats": "csv,tsv,json",
            "pseudonymize": None,
            "expect_exit": "zero",
            "allow_missing_inputs": False,
        },
        {
            "case_id": "cohort_empty",
            "kind": "cohort",
            "group": "cohort",
            "inputs": [],
            "empty_input_dir": True,
            "summary_formats": "csv,tsv,json",
            "pseudonymize": None,
            "expect_exit": "zero",
            "allow_missing_inputs": True,
        },
    ]


def check_matrix(cases: list[dict[str, Any]], probes: list[dict[str, Any]]) -> dict[str, Any]:
    """Compare the derived matrix against the counts the gate page records.

    Args:
        cases: All pipeline cases.
        probes: The probe cases.

    Returns:
        dict[str, Any]: ``counts``, ``documented`` and ``mismatches``. An empty
        ``mismatches`` means this run's matrix is the one runs 1-3 measured.
    """
    by_assembly: dict[str, int] = {}
    by_group: dict[str, int] = {}
    for case in cases:
        by_assembly[case["assembly"]] = by_assembly.get(case["assembly"], 0) + 1
        by_group[case["group"]] = by_group.get(case["group"], 0) + 1

    assembly_counts = dict(sorted(by_assembly.items()))
    totals: dict[str, int] = {
        "base": by_group.get("base", 0),
        "nonfast": by_group.get("nonfast", 0),
        "advntr": by_group.get("advntr", 0),
        "total": len(cases),
        "probes": len(probes),
    }
    counts: dict[str, Any] = {"by_assembly": assembly_counts, **totals}
    documented: dict[str, Any] = {"by_assembly": DOCUMENTED_ASSEMBLY_COUNTS, **DOCUMENTED_TOTALS}

    mismatches: list[str] = []
    for key, expected in DOCUMENTED_TOTALS.items():
        if totals[key] != expected:
            mismatches.append(f"{key}: derived {totals[key]}, page records {expected}")
    for assembly, expected in DOCUMENTED_ASSEMBLY_COUNTS.items():
        actual = assembly_counts.get(assembly, 0)
        if actual != expected:
            mismatches.append(f"assembly {assembly}: derived {actual}, page records {expected}")
    mismatches.extend(
        f"assembly {assembly}: derived {count}, page records none"
        for assembly, count in assembly_counts.items()
        if assembly not in DOCUMENTED_ASSEMBLY_COUNTS
    )

    return {"counts": counts, "documented": documented, "mismatches": mismatches}


def build_matrix(
    data_dir: Path,
    *,
    non_fast_ids: tuple[str, ...] = NON_FAST_CASE_IDS,
    advntr_ids: tuple[str, ...] = ADVNTR_CASE_IDS,
    advntr_max_coverage: int = 300,
    case_filter: list[str] | None = None,
    include_probes: bool = True,
    include_cohort: bool = True,
) -> dict[str, Any]:
    """Derive the whole matrix, log what was derived, and self-check it against the page.

    Args:
        data_dir: The ``tests/data`` directory.
        non_fast_ids: Override for :data:`NON_FAST_CASE_IDS`.
        advntr_ids: Override for :data:`ADVNTR_CASE_IDS`.
        advntr_max_coverage: The ``--advntr-max-coverage`` value.
        case_filter: If given, keep only cases whose id contains one of these substrings.
            For smoke tests; the check against the page's counts is skipped when it is used
            and the matrix records that it was.
        include_probes: Whether to build the deliberate-mismatch probes.
        include_cohort: Whether to build the cohort-mode cases.

    Returns:
        dict[str, Any]: The matrix, ready to be written to ``matrix.json``.
    """
    base_cases, log = derive_base_cases(data_dir)
    cases, policy_log = apply_policies(
        base_cases,
        non_fast_ids=non_fast_ids,
        advntr_ids=advntr_ids,
        advntr_max_coverage=advntr_max_coverage,
    )
    log.extend(policy_log)
    probes = build_probes(base_cases) if include_probes else []

    if case_filter:
        kept = [case for case in cases if any(needle in case["case_id"] for needle in case_filter)]
        kept_probes = [probe for probe in probes if any(needle in probe["case_id"] for needle in case_filter)]
        log.append(
            f"case filter {case_filter!r}: {len(kept)}/{len(cases)} cases, {len(kept_probes)}/{len(probes)} probes"
        )
        cases, probes = kept, kept_probes

    cases.sort(key=lambda case: (case["group"], case["case_id"]))
    probes.sort(key=lambda case: case["case_id"])

    cohort_cases = build_cohort_cases([case["case_id"] for case in cases]) if include_cohort else []
    check = check_matrix(cases, probes)
    check["skipped"] = bool(case_filter)

    for line in log:
        logger.info(f"matrix: {line}")
    logger.info(f"matrix: counts {check['counts']}")
    if check["skipped"]:
        logger.warning("matrix: a case filter is in force, so the check against the gate page's counts is advisory")
    for mismatch in check["mismatches"]:
        logger.warning(f"matrix: DIFFERS FROM THE GATE PAGE - {mismatch}")

    return {
        "data_dir": str(data_dir),
        "derivation_log": log,
        "policies": {
            "non_fast_case_ids": list(non_fast_ids),
            "advntr_case_ids": list(advntr_ids),
            "advntr_max_coverage": advntr_max_coverage,
            "probe_specs": [list(spec) for spec in PROBE_SPECS] if include_probes else [],
        },
        "check": check,
        "cases": cases,
        "probes": probes,
        "cohort_cases": cohort_cases,
    }
