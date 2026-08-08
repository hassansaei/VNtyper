"""Focused CRAM-case policy and fixture-path derivation for the golden cohort."""

from __future__ import annotations

import logging
from collections.abc import Callable
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

CRAM_SCAN_MODES: tuple[str, ...] = ("indexed", "stream")


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
            cases.append(case)

    log.append(f"CRAM repeats: {len(cases)}/{len(cram_ids) * len(CRAM_SCAN_MODES)} from {cram_root}")
    if missing:
        log.append(f"CRAM fixtures missing for: {', '.join(missing)}")
    return cases, log
