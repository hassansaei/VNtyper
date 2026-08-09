"""Focused CRAM-case policy and fixture-path derivation for the golden cohort."""

from __future__ import annotations

import logging
from collections.abc import Callable
from pathlib import Path
from typing import Any

from golden_cohort.admissibility import PIPELINE_REQUIRED_ARTIFACTS

logger = logging.getLogger(__name__)

CRAM_SCAN_MODES: tuple[str, ...] = ("indexed", "stream")

CRAM_READ_SET_EXPECTATIONS: dict[str, dict[str, Any]] = {
    "7a61_hg38_ensembl_bwa": {
        "placed_unmapped_guard_count": 11_571,
        "raw_indexed_read_set": {
            "count": 2_690,
            "sorted_read_name_sha256": "c64be739cd6344b8b62004fc9ea568779b3cc06ff1d472ac0e5d97c343130d7d",
        },
        "stream_read_set": {
            "count": 634_261,
            "sorted_read_name_sha256": "b7f75d19497698f12d6dbbc38afc12702b2d262670a4c893b39f95967ebf7b7b",
        },
    },
    "b178_hg19_subset": {
        "placed_unmapped_guard_count": 329,
        "raw_indexed_read_set": {
            "count": 4_478,
            "sorted_read_name_sha256": "dad9a81a4e8cf30d1d938717459614f7d8ac6decb84978a5bc23c090b4d90a8b",
        },
        "stream_read_set": {
            "count": 4_807,
            "sorted_read_name_sha256": "d3aa88fe91c8964b2f9a1b053a672f2bc3d1896b71de986f5cde02999d552591",
        },
    },
}

INDEXED_SAFE_READ_SET: dict[str, int | str] = {
    "count": 20,
    "sorted_read_name_sha256": "16a0efa7785630c3d80716d9a386ddaa24f4933b5671f4ecd221b42a8dffe740",
}

# A no-variant pipeline success does not write the pre-filter variant frame. The purpose
# cases still require the final Kestrel table, coverage, pipeline summary and HTML report.
INDEXED_SAFE_REQUIRED_ARTIFACTS: tuple[str, ...] = tuple(
    artifact for artifact in PIPELINE_REQUIRED_ARTIFACTS if artifact != "kestrel/kestrel_pre_result.tsv"
)


def _indexed_safe_cases(fixture: Path) -> list[dict[str, Any]]:
    """Build the nonempty purpose cases that authorize both extraction strategies."""
    expectation: dict[str, Any] = {
        "indexed_authorized": True,
        "raw_indexed_read_set": dict(INDEXED_SAFE_READ_SET),
        "stream_read_set": dict(INDEXED_SAFE_READ_SET),
    }
    return [
        {
            "case_id": f"indexed_safe_{scan}_cram",
            "kind": "pipeline",
            "group": "cram",
            "alignment_kind": "cram",
            "sample": "example_indexed_safe",
            "assembly": "hg19",
            "cram": str(fixture),
            "fast_mode": False,
            "advntr": False,
            "unmapped_scan": scan,
            "repeat_of": "indexed_safe",
            "expect_exit": "zero",
            "required_artifacts": list(INDEXED_SAFE_REQUIRED_ARTIFACTS),
            "cram_evidence_expectation": expectation,
        }
        for scan in CRAM_SCAN_MODES
    ]


def cram_fixture_for(case: dict[str, Any], data_dir: Path, cram_root: Path) -> Path:
    """Return the mirrored CRAM path for a derived base BAM case.

    Args:
        case: A derived base case with a BAM under ``data_dir``.
        data_dir: The test-data directory.
        cram_root: The root written by the fixture deriver.

    Returns:
        The corresponding CRAM path.

    Raises:
        ValueError: If the case BAM is not under ``data_dir``.
    """
    source = Path(case["bam"])
    try:
        relative = source.relative_to(data_dir.resolve())
    except ValueError:
        msg = (
            f"Cannot derive a CRAM fixture path for {case['case_id']}: its BAM {source} is not under "
            f"the data directory {data_dir}. The fixture tree mirrors the data directory, so a BAM "
            "from outside it has no mirrored position."
        )
        logger.error(msg)
        raise ValueError(msg) from None
    return (cram_root / relative).with_suffix(".cram")


def build_cram_cases(
    base_cases: list[dict[str, Any]],
    *,
    cram_ids: tuple[str, ...],
    data_dir: Path,
    cram_root: Path,
    resolve: Callable[[tuple[str, ...], dict[str, dict[str, Any]], str], list[dict[str, Any]]],
) -> tuple[list[dict[str, Any]], list[str]]:
    """Build non-fast CRAM repeats for each lossless scan strategy.

    Args:
        base_cases: Derived BAM-by-assembly cases.
        cram_ids: Policy-selected base case ids.
        data_dir: Test-data root used to mirror source paths.
        cram_root: Root containing derived CRAM fixtures.
        resolve: Policy resolver supplied by the enclosing matrix module.

    Returns:
        CRAM cases and derivation-log lines. Missing fixtures remain visible as matrix
        drift rather than becoming a reduced successful matrix.
    """
    by_id = {case["case_id"]: case for case in base_cases}
    cases: list[dict[str, Any]] = []
    log: list[str] = []
    missing: list[str] = []

    for base in resolve(cram_ids, by_id, "CRAM"):
        fixture = cram_fixture_for(base, data_dir, cram_root)
        if not fixture.is_file():
            missing.append(base["case_id"])
            log.append(f"skipped (no derived CRAM fixture at {fixture}): {base['case_id']}")
            logger.error(
                f"matrix: no CRAM fixture for {base['case_id']} at {fixture}. Run `make cram-fixtures` to "
                "derive it. This run's CRAM group is short, which the matrix check reports as drift."
            )
            continue
        for scan in CRAM_SCAN_MODES:
            case = dict(base)
            case.pop("bam", None)
            evidence_expectation = CRAM_READ_SET_EXPECTATIONS.get(base["case_id"])
            case.update(
                {
                    "case_id": f"{base['sample'].removeprefix('example_')}_{base['assembly']}_{scan}_cram",
                    "group": "cram",
                    "alignment_kind": "cram",
                    "cram": str(fixture),
                    "source_bam": base["bam"],
                    "fast_mode": False,
                    "unmapped_scan": scan,
                    "repeat_of": base["case_id"],
                }
            )
            if evidence_expectation is not None:
                case["side_expectations"] = {
                    "before": {
                        "expect_exit": "zero",
                        "required_artifacts": list(base["required_artifacts"]),
                    },
                    "after": {
                        "expect_exit": "nonzero",
                        "required_artifacts": [],
                        "cram_evidence_expectation": evidence_expectation,
                    },
                }
            cases.append(case)

    if cram_ids:
        indexed_safe_fixture = cram_root / "indexed-safe" / "indexed-safe.cram"
        if indexed_safe_fixture.is_file():
            cases.extend(_indexed_safe_cases(indexed_safe_fixture))
        else:
            missing.append("indexed_safe")
            log.append(f"skipped (no indexed-safe CRAM fixture at {indexed_safe_fixture}): indexed_safe")
            logger.error(
                f"matrix: no indexed-safe CRAM fixture at {indexed_safe_fixture}. Run `make cram-fixtures` to "
                "derive it. This run's CRAM group is short, which the matrix check reports as drift."
            )

    declared_sources = len(cram_ids) + (1 if cram_ids else 0)
    log.append(f"CRAM repeats: {len(cases)}/{declared_sources * len(CRAM_SCAN_MODES)} from {cram_root}")
    if missing:
        log.append(f"CRAM fixtures missing for: {', '.join(missing)}")
    return cases, log
