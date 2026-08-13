"""Performance budgets for the nomenclature path (spec §3.4, §3.5).

These live in the unit tier, marked ``unit``, so ``make check-all`` runs them: a
budget that does not gate is not a budget.

The thresholds are deliberately the spec's, not the measured figures. Measured on
this machine the per-call cost is roughly an order of magnitude under budget, so
the headroom absorbs a loaded CI box without the assertion becoming decorative.

Research use only.
"""

from __future__ import annotations

import importlib
import time
from pathlib import Path

import pysam
import pytest

from vntyper.scripts.nomenclature import from_advntr, from_kestrel
from vntyper.scripts.nomenclature_bam import BamRescuer, is_candidate

pytestmark = pytest.mark.unit

#: Spec §3.4.
TABLE_BUILD_BUDGET_MS = 20.0
PER_CALL_BUDGET_US = 50.0
#: Spec §3.5.
PER_LOCUS_BUDGET_MS = 15.0


def _write_bam(path: Path, reads: int = 200) -> Path:
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "K-J", "LN": 120}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as handle:
        for index in range(reads):
            record = pysam.AlignedSegment(handle.header)
            record.query_name = f"r{index}"
            record.reference_id = 0
            record.reference_start = 10
            record.mapping_quality = 255
            record.cigarstring = "20=1X1I20=" if index % 4 == 0 else "20=1X20="
            record.query_sequence = "A" * (42 if index % 4 == 0 else 41)
            handle.write(record)
    pysam.index(str(path))  # type: ignore[attr-defined]
    return path


def test_the_lookup_tables_build_within_budget() -> None:
    """Everything is built once at import; a slow import is paid by every process."""
    import vntyper.scripts.nomenclature as module

    start = time.perf_counter()
    importlib.reload(module)
    elapsed_ms = (time.perf_counter() - start) * 1000
    assert elapsed_ms < TABLE_BUILD_BUDGET_MS, f"table build took {elapsed_ms:.1f} ms"


def test_from_kestrel_stays_within_the_per_call_budget() -> None:
    """It is called from a dataframe apply, so per-call cost is the whole cost."""
    from_kestrel("S-C", 67, "G", "GG")  # warm

    iterations = 2000
    start = time.perf_counter()
    for _ in range(iterations):
        from_kestrel("S-C", 67, "G", "GG")
    per_call_us = (time.perf_counter() - start) / iterations * 1_000_000
    assert per_call_us < PER_CALL_BUDGET_US, f"{per_call_us:.1f} us/call"


def test_from_advntr_stays_within_the_per_call_budget() -> None:
    """The state parser is compiled once at import, never in the hot path."""
    from_advntr("D27_2&I27_2_A_LEN2")  # warm

    iterations = 2000
    start = time.perf_counter()
    for _ in range(iterations):
        from_advntr("D27_2&I27_2_A_LEN2")
    per_call_us = (time.perf_counter() - start) / iterations * 1_000_000
    assert per_call_us < PER_CALL_BUDGET_US, f"{per_call_us:.1f} us/call"


def test_a_candidate_locus_is_rescued_within_budget(tmp_path: Path) -> None:
    bam = _write_bam(tmp_path / "perf.bam")
    rescuer = BamRescuer(bam)
    rescuer.rescue("K-J", 30)  # warm, and opens the handle once

    iterations = 20
    start = time.perf_counter()
    for _ in range(iterations):
        rescuer.rescue("K-J", 30)
    per_locus_ms = (time.perf_counter() - start) / iterations * 1000
    assert per_locus_ms < PER_LOCUS_BUDGET_MS, f"{per_locus_ms:.2f} ms/locus"


def test_the_rescue_path_does_not_run_when_no_candidate_exists(tmp_path: Path) -> None:
    """The spec's other half of the budget: with no candidate the cost is zero.

    Asserted structurally rather than by timing -- the handle is never opened and
    no region is ever fetched.
    """
    from vntyper.scripts.nomenclature import Nomenclature

    confident = Nomenclature(
        name="59dupC",
        event="duplication",
        unit="X",
        tier="A",
        flags=(),
        ambiguity=(53, 59),
        repeat_form="53C[7]>53C[8]",
        net_length=1,
        source="reconciled",
    )
    rescuer = BamRescuer(_write_bam(tmp_path / "unused.bam"))
    if is_candidate(confident):
        rescuer.rescue("K-J", 30)

    assert rescuer.opens == 0
    assert rescuer.fetches == 0
    assert rescuer.opened is False
