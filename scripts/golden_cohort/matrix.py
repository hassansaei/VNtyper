"""Derive the golden-cohort case matrix from what is actually in ``tests/data``.

The gate page describes the matrix in prose: "the 7 multi-reference samples at all six
assemblies plus their original hg19 subsets, and the hg38 regression guard
``example_40cf``", then five repeats without ``--fast-mode`` and three ``--extra-modules
advntr`` runs. A hardcoded list of ids reproduces that only until ``tests/data``
changes, and then reproduces it wrongly and silently. So the 50 BAM-by-assembly cases are
**derived** by walking the data directory, and the four selections that are policy rather
than data - which cases repeat without fast mode, which run adVNTR, which repeat from a
derived CRAM, which are the deliberate-mismatch probes - are declared here by case id and
**resolved** against the derived set. A policy naming a case the data does not contain is
an error, not a silent drop.

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

The CRAM group
--------------
VNtyper accepts CRAM and, up to and including run 5, no gate run had ever fed it one - so
the CRAM branch of ``process_bam_to_fastq`` and the write race ``175011e`` fixed in
``build_cram_unmapped_filter_command`` were attested by unit tests and one hand-run
equivalence comparison, never by this gate. ``make cram-fixtures``
(``scripts/make_cram_fixtures.py``) now derives a verified CRAM beside every cohort BAM
under ``tests/data/cram/``, mirroring the source layout with ``.bam`` -> ``.cram``, and
:data:`CRAM_CASE_IDS` declares which of the derived base cases repeat from that CRAM.

Two properties of these cases are load-bearing rather than incidental:

* they run **without** ``--fast-mode``. ``fast_mode=True`` skips the unmapped-read
  extraction entirely (``fastq_bam_processing.process_bam_to_fastq``, ``if not
  fast_mode:``), and the CRAM-specific ``build_cram_unmapped_filter_command`` branch lives
  inside it. A fast-mode CRAM case would exercise the slice and the FASTQ conversion and
  none of the code the fixtures exist for.
* a declared CRAM case whose fixture has not been derived is **skipped and logged**, never
  silently dropped from the contract. The count then comes out short, which is an ordinary
  :func:`check_matrix` mismatch and refuses the run unless ``--allow-matrix-drift`` is
  passed. There is deliberately no "0 or 6 CRAM cases are both fine" rule: a run without
  them is a reduced run and must not earn an attestation-grade verdict.

What "derived" does and does not mean
-------------------------------------
Only the **base cases** are derived. The five non-fast ids, the three adVNTR ids, the two
cohort CRAM ids, the indexed-safe purpose CRAM and the three probes are hardcoded policy,
with the base-case selections resolved against the derived set. The
CRAM *fixture paths* are derived from the base case's BAM path, so they cannot drift from
what ``make_cram_fixtures.py`` wrote, but which cases are chosen is policy like the rest.
Anything that describes this matrix as "derived" without that qualification is overstating
it, and the gate page has done exactly that.

Drift is fatal by default
-------------------------
:func:`check_matrix` compares the derivation against the per-group contract the gate page
records - 50 base, 5 non-fast, 3 adVNTR, 6 CRAM and 3 probes. (It was 50/5/3 plus 3 probes
for runs 1-5, which is the matrix every result table on that page was measured over; run 6
took an earlier, smaller CRAM group.) That check used to be advisory in every
direction: ``build_matrix`` logged the
deviations as warnings, ``cmd_matrix`` returned 0 regardless, and the comparison's verdict
ignored them - so a silently reduced run earned the same ``IDENTICAL`` as a full one, and a
``--case`` filter matching nothing produced a zero-case matrix that every ``all()`` in the
harness then agreed was verified. ``strict=True`` now refuses to build a drifted matrix at
all, and a zero-case matrix is refused whether strict or not. A deliberately reduced run is
still available through ``--case`` or ``--allow-matrix-drift``, but it is recorded as not
attestation-grade and the comparison gives it a different verdict word.

Attributes:
    ASSEMBLIES: The six assemblies the cohort is provided at, in the page's order.
    NON_FAST_CASE_IDS: Which derived cases repeat without ``--fast-mode``.
    ADVNTR_CASE_IDS: Which derived cases repeat with ``--extra-modules advntr``.
    CRAM_CASE_IDS: Which derived cases repeat from their derived CRAM fixture.
    CRAM_FIXTURE_DIRNAME: The fixture root, relative to the data directory.
    PROBE_SPECS: The three deliberate-mismatch probes.
    DOCUMENTED_ASSEMBLY_COUNTS: The per-assembly case counts, used as a self-check.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any

from golden_cohort.admissibility import PIPELINE_REQUIRED_ARTIFACTS
from golden_cohort.case_expectations import (
    declare_mixed_layout_outcome,
    materialize_side_expectation,
    without_side_expectations,
)
from golden_cohort.cohort_cases import build_cohort_cases
from golden_cohort.cram_cases import build_cram_cases

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
#: ``dfc3_hg19_advntr``, all at ``--advntr-max-coverage 300``. ``b178`` uses the
#: measured clean remapped BAM so candidate cohort runs retain a real adVNTR producer.
ADVNTR_CASE_IDS: tuple[str, ...] = (
    "a5c1_hg19_subset",
    "b178_hg19_bwa",
    "dfc3_hg19_subset",
)

#: Which base cases repeat from their derived CRAM fixture. Policy, like the two selections
#: above, and chosen to cover both a call and the write race:
#:
#: * ``b178_hg19_subset`` - a known historical positive (``D-C`` insertion,
#:   ``High_Precision*``, not flagged, per the gate page's run-1 table), 34,214 records
#:   and 4,478 flag-12 pairs in the fixture manifest. Corrected flag-4 extraction records
#:   4,807 reads and exposes the 329 reads an indexed ``'*'`` fetch would lose.
#: * ``7a61_hg38_ensembl_bwa`` - 985,731 records and 622,690 flag-12 pairs, one of the
#:   heaviest unmapped loads in the cohort, and so the most exposed to the write race
#:   ``175011e`` fixed (measured there at 199,797 of 200,000 unmapped reads present when the
#:   shell returned). Corrected flag-4 extraction records 634,261 reads.
#:
#: It is **not** the single heaviest: ``7a61_hg19_subset`` carries 958,804 unmapped pairs,
#: and the six remapped ``7a61`` cases tie at 623,792 / 622,690. This pair is kept because it
#: is the pair already proven end to end, and because the two ids sit at the two different
#: layouts the fixture tree mirrors - one top-level subset BAM and one
#: ``remapped/<aligner>/<assembly>/`` BAM - so the path derivation below is exercised on
#: both shapes rather than only one. Counts are from ``tests/data/cram/manifest.json``.
CRAM_CASE_IDS: tuple[str, ...] = (
    "b178_hg19_subset",
    "7a61_hg38_ensembl_bwa",
)

#: Where ``scripts/make_cram_fixtures.py`` writes, relative to the data directory. It
#: mirrors the source layout underneath, so a fixture path is derived from its BAM path
#: rather than declared.
CRAM_FIXTURE_DIRNAME = "cram"

#: ``(probe_id, base_case_id, declared_assembly, expectation)``. Two deliberate
#: mismatches, which the page records as exit 1 on both sides with only the failure point
#: moving, and one naming probe (an NCBI-named BAM declared ``hg38``) which exits 0.
PROBE_SPECS: tuple[tuple[str, str, str, str], ...] = (
    ("probe_mismatch_hg19_as_hg38", "dfc3_hg19_subset", "hg38", "nonzero"),
    ("probe_mismatch_hg38_as_hg19", "40cf_hg38_subset", "hg19", "nonzero"),
    ("probe_naming_ncbi_as_hg38", "dfc3_GRCh38_bwa", "hg38", "zero"),
)

#: The per-assembly case counts, used as a self-check on the derivation.
#:
#: Runs 1-5 measured 20 hg19 / 9 hg38 / 8 GRCh38 / 7 GRCh37 / 7 hg19_ensembl /
#: 7 hg38_ensembl, and run 2's after-side assembly-guard verdict counts are that same
#: distribution. Each selected CRAM fixture now runs through both lossless scan strategies,
#: moving hg19 to 24 and hg38_ensembl to 9. The page's ``x / 58`` result tables therefore
#: describe the pre-CRAM matrix and are not restated here.
DOCUMENTED_ASSEMBLY_COUNTS: dict[str, int] = {
    "hg19": 24,
    "hg38": 9,
    "GRCh38": 8,
    "GRCh37": 7,
    "hg19_ensembl": 7,
    "hg38_ensembl": 9,
}

#: The page's own totals, checked the same way. ``total`` was 58 for runs 1-5; six CRAM
#: cases cover indexed and stream extraction and take it to 64.
DOCUMENTED_TOTALS: dict[str, int] = {"base": 50, "nonfast": 5, "advntr": 3, "cram": 6, "total": 64, "probes": 3}


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
            declare_mixed_layout_outcome(
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
                    "required_artifacts": list(PIPELINE_REQUIRED_ARTIFACTS),
                }
            )
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
                declare_mixed_layout_outcome(
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
                        "required_artifacts": list(PIPELINE_REQUIRED_ARTIFACTS),
                    }
                )
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
        case = without_side_expectations(base)
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

    A probe expected to exit nonzero declares **no** required artefacts: the whole point of
    the two mismatch probes is that the run refuses, and run 4 measured that they produce
    none of the five pipeline artefacts on either side. The naming probe is expected to
    exit zero and therefore carries the full requirement, which it met on both sides.

    Returns:
        list[dict]: The probe cases. Probes are run and compared but are not part of the
        58-case matrix, exactly as the page counts them.
    """
    by_id = {case["case_id"]: case for case in base_cases}
    probes: list[dict[str, Any]] = []
    for probe_id, base_id, declared, expectation in PROBE_SPECS:
        base = _resolve((base_id,), by_id, f"probe {probe_id}")[0]
        probe = without_side_expectations(base)
        probe.update(
            {
                "case_id": probe_id,
                "group": "probe",
                "assembly": declared,
                "declared_assembly": declared,
                "true_assembly": base["assembly"],
                "expect_exit": expectation,
                "repeat_of": base["case_id"],
                "required_artifacts": list(PIPELINE_REQUIRED_ARTIFACTS) if expectation == "zero" else [],
            }
        )
        probes.append(probe)
    return probes


def check_matrix(cases: list[dict[str, Any]], probes: list[dict[str, Any]]) -> dict[str, Any]:
    """Compare the derived matrix against the counts the gate page records.

    Args:
        cases: All pipeline cases.
        probes: The probe cases.

    Returns:
        dict[str, Any]: ``counts``, ``documented``, ``mismatches`` and
        ``attestation_grade``. An empty ``mismatches`` means this run's matrix is the one
        the gate page's contract describes; ``attestation_grade`` additionally requires that
        no case filter narrowed it, and is what :func:`golden_cohort.compare._verdict` reads
        to decide whether a clean result may be called ``IDENTICAL`` at all. It is
        deliberately *not* the matrix runs 1-5 measured: those predate the CRAM group, so a
        clean run of this matrix covers strictly more than they did.
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
        # Counted like the rest rather than exempted: a run that derived no CRAM fixtures
        # is a reduced run, and reporting 0 here is what makes it say so.
        "cram": by_group.get("cram", 0),
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

    return {
        "counts": counts,
        "documented": documented,
        "mismatches": mismatches,
        "skipped": False,
        "attestation_grade": not mismatches,
    }


def build_matrix(
    data_dir: Path,
    *,
    non_fast_ids: tuple[str, ...] = NON_FAST_CASE_IDS,
    advntr_ids: tuple[str, ...] = ADVNTR_CASE_IDS,
    cram_ids: tuple[str, ...] = CRAM_CASE_IDS,
    cram_root: Path | None = None,
    advntr_max_coverage: int = 300,
    case_filter: list[str] | None = None,
    include_probes: bool = True,
    include_cohort: bool = True,
    strict: bool = True,
) -> dict[str, Any]:
    """Derive the whole matrix, log what was derived, and self-check it against the page.

    Args:
        data_dir: The ``tests/data`` directory.
        non_fast_ids: Override for :data:`NON_FAST_CASE_IDS`.
        advntr_ids: Override for :data:`ADVNTR_CASE_IDS`.
        cram_ids: Override for :data:`CRAM_CASE_IDS`.
        cram_root: Where ``make cram-fixtures`` wrote. Defaults to
            ``data_dir / CRAM_FIXTURE_DIRNAME``.
        advntr_max_coverage: The ``--advntr-max-coverage`` value.
        case_filter: If given, keep only cases whose id contains one of these substrings.
            For smoke tests; the check against the page's counts becomes advisory when it
            is used, the matrix records that it was, and the result is marked as not
            attestation-grade so the comparison cannot call it ``IDENTICAL``.
        include_probes: Whether to build the deliberate-mismatch probes.
        include_cohort: Whether to build the cohort-mode cases.
        strict: Refuse to return an unfiltered matrix that deviates from the documented
            contract. Set False only for a deliberately reduced run, which is then also
            marked as not attestation-grade.

    Returns:
        dict[str, Any]: The matrix, ready to be written to ``matrix.json``.

    Raises:
        ValueError: If the matrix has no pipeline cases at all, or if ``strict`` and no
            case filter is in force and the derivation deviates from the documented
            counts. Both used to be warnings. A zero-case matrix is the more insidious of
            the two: ``all()`` over an empty mapping is True, so a side that ran nothing
            reported ``launch_verified`` and compared clean against another side that also
            ran nothing.
    """
    base_cases, log = derive_base_cases(data_dir)
    cases, policy_log = apply_policies(
        base_cases,
        non_fast_ids=non_fast_ids,
        advntr_ids=advntr_ids,
        advntr_max_coverage=advntr_max_coverage,
    )
    log.extend(policy_log)

    resolved_cram_root = cram_root if cram_root is not None else data_dir / CRAM_FIXTURE_DIRNAME
    cram_cases, cram_log = build_cram_cases(
        base_cases,
        cram_ids=cram_ids,
        data_dir=data_dir,
        cram_root=resolved_cram_root,
        resolve=_resolve,
    )
    cases.extend(cram_cases)
    log.extend(cram_log)

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

    candidate_output_ids = [
        case["case_id"] for case in cases if materialize_side_expectation(case, "after").get("expect_exit") == "zero"
    ]
    cohort_cases = build_cohort_cases(candidate_output_ids) if include_cohort else []
    check = check_matrix(cases, probes)
    check["skipped"] = bool(case_filter)
    check["strict"] = bool(strict)
    check["attestation_grade"] = not check["mismatches"] and not check["skipped"]

    for line in log:
        logger.info(f"matrix: {line}")
    logger.info(f"matrix: counts {check['counts']}")

    if not cases:
        msg = (
            f"The matrix has no pipeline cases (data_dir={data_dir}, case filter={case_filter!r}). "
            "Refusing to build it: a zero-case side runs nothing, and every `all()` in this harness "
            "agrees vacuously that nothing was verified, so the comparison would report IDENTICAL over "
            "two runs that never happened."
        )
        logger.error(msg)
        raise ValueError(msg)

    if check["skipped"]:
        logger.warning(
            "matrix: a case filter is in force, so this run is NOT attestation-grade and its comparison "
            "verdict will say so"
        )
    for mismatch in check["mismatches"]:
        logger.warning(f"matrix: DIFFERS FROM THE GATE PAGE - {mismatch}")
    if check["mismatches"] and strict and not check["skipped"]:
        msg = (
            f"The derived matrix deviates from the contract the gate page records: "
            f"{'; '.join(check['mismatches'])}. Refusing to launch an attestation run over a matrix that is "
            "not the one the page describes - a reduced run earns the same IDENTICAL verdict as a full one, "
            "which is how a shrinking gate stays invisible. Either fix tests/data or the policy, or pass "
            "--allow-matrix-drift to run it knowingly as a non-attestation run."
        )
        logger.error(msg)
        raise ValueError(msg)

    return {
        "data_dir": str(data_dir),
        "derivation_log": log,
        "policies": {
            "non_fast_case_ids": list(non_fast_ids),
            "advntr_case_ids": list(advntr_ids),
            "cram_case_ids": list(cram_ids),
            "cram_root": str(resolved_cram_root),
            "advntr_max_coverage": advntr_max_coverage,
            "probe_specs": [list(spec) for spec in PROBE_SPECS] if include_probes else [],
        },
        "check": check,
        "cases": cases,
        "probes": probes,
        "cohort_cases": cohort_cases,
    }
