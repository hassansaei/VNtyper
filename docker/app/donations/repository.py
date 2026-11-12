"""Storage and k-suppressed aggregates for anonymous research donations."""

from __future__ import annotations

import logging
import sqlite3
import uuid
from collections import Counter
from datetime import datetime, timezone
from typing import Any

from ..config import settings
from .models import KITS, SCHEMA, TRUTHS, Donation

logger = logging.getLogger(__name__)

K_DEFAULT = 5
# Cells that count confirmed positives or tool calls describe a rare disease.
K_SENSITIVE = 10

# Whitelist columns, in Donation field order; `schema` is implied by the table name.
COLUMNS = tuple(name for name in Donation.model_json_schema()["properties"] if name != "schema")

# `id` is random and never leaves the server.
CREATE_SQL = [
    """CREATE TABLE IF NOT EXISTS donations_v2 (
        id TEXT PRIMARY KEY,
        received_quarter TEXT NOT NULL,
        vntyper_version TEXT NOT NULL,
        coverage_qc TEXT NOT NULL,
        flank_depth_bin TEXT NOT NULL,
        array_depth_bin TEXT NOT NULL,
        kestrel_call BOOLEAN NOT NULL,
        kestrel_flagged BOOLEAN NOT NULL,
        kestrel_confidence TEXT NOT NULL,
        depth_score_bin TEXT,
        variant_class TEXT NOT NULL,
        kit TEXT NOT NULL,
        platform_class TEXT NOT NULL,
        truth TEXT NOT NULL,
        confirmation_class TEXT NOT NULL,
        referral TEXT
    )""",
]


def received_quarter(now: datetime | None = None) -> str:
    """Return the current UTC quarter as ``YYYY-Qn``, the finest time a row keeps."""
    now = now or datetime.now(timezone.utc)
    return f"{now.year}-Q{(now.month - 1) // 3 + 1}"


def suppress(counts: dict[str, int], k: dict[str, int]) -> dict[str, int | None]:
    """Hide cells below their k, then hide more until the hidden set cannot be recovered.

    Primary suppression hides each non-zero cell under its threshold (a zero names nobody).
    Complementary suppression then hides the next smallest non-zero cells until at least
    two are hidden and they sum to at least the largest threshold among them, so
    ``total - visible`` never pins down a single small cell.
    """
    hidden = {c for c, n in counts.items() if 0 < n < k.get(c, K_DEFAULT)}

    def recoverable() -> bool:
        need = max(k.get(c, K_DEFAULT) for c in hidden)
        return len(hidden) < 2 or sum(counts[c] for c in hidden) < need

    while hidden and recoverable():
        spare = [c for c in counts if c not in hidden and counts[c] > 0]
        if not spare:
            break
        hidden.add(min(spare, key=lambda c: counts[c]))
    return {c: None if c in hidden else n for c, n in counts.items()}


class DonationRepository:
    def __init__(self, db_url: str | None = None) -> None:
        url = db_url or settings.POSTGRES_URL
        if not url:
            raise RuntimeError("Donations are enabled but POSTGRES_URL is not set; refusing to use a throwaway store.")
        self.db_url: str = url
        self._is_postgres = url.startswith(("postgresql://", "postgres://"))
        self._param = "%s" if self._is_postgres else "?"
        self._execute(CREATE_SQL)
        self._drop_v1_table_if_empty()

    def _drop_v1_table_if_empty(self) -> None:
        """Drop the v1 table (run IDs, digests, exact timestamps) only when it holds nothing."""
        exists_sql = (
            "SELECT COUNT(*) FROM information_schema.tables WHERE table_name = 'donations'"
            if self._is_postgres
            else "SELECT COUNT(*) FROM sqlite_master WHERE type = 'table' AND name = 'donations'"
        )
        if not self._execute([exists_sql])[0][0]:
            return
        rows = self._execute(["SELECT COUNT(*) FROM donations"])[0][0]
        if rows:
            logger.error(
                "Legacy table 'donations' holds %d pseudonymous v1 rows and was left in place. "
                "Export or delete them manually, then drop the table.",
                rows,
            )
        else:
            self._execute(["DROP TABLE donations"])

    def _connect(self) -> Any:
        if self._is_postgres:
            import psycopg

            return psycopg.connect(self.db_url)
        return sqlite3.connect(self.db_url)  # a file path; tests use this

    def _execute(self, statements: list[str], params: tuple[Any, ...] = ()) -> list[Any]:
        conn = self._connect()
        try:
            cur = conn.cursor()
            for sql in statements:
                cur.execute(sql, params) if params else cur.execute(sql)
            rows = cur.fetchall() if cur.description else []
            conn.commit()
            return list(rows)
        finally:
            conn.close()

    def save(self, donation: Donation) -> None:
        values = donation.model_dump(exclude={"schema_"})
        cols = ("id", "received_quarter", *COLUMNS)
        sql = f"INSERT INTO donations_v2 ({', '.join(cols)}) VALUES ({', '.join([self._param] * len(cols))})"
        self._execute([sql], (str(uuid.uuid4()), received_quarter(), *(values[c] for c in COLUMNS)))

    def get_aggregates(self) -> dict[str, Any]:
        # Closed quarters only: a live count would reveal each new donation by differencing.
        rows = self._execute(
            [
                f"SELECT kit, truth, kestrel_call, COUNT(*) FROM donations_v2 WHERE received_quarter < {self._param} "
                "GROUP BY kit, truth, kestrel_call"
            ],
            (received_quarter(),),
        )
        by_kit: Counter[str] = Counter()
        by_truth: Counter[str] = Counter()
        by_call: Counter[str] = Counter()
        for kit, truth, call, n in rows:
            by_kit[kit] += n
            by_truth[truth] += n
            by_call["true" if call else "false"] += n
        total = sum(by_kit.values())
        tables = {
            "by_truth": suppress({t: by_truth[t] for t in TRUTHS}, {"confirmed_positive": K_SENSITIVE}),
            "by_kestrel_call": suppress({c: by_call[c] for c in ("true", "false")}, {"true": K_SENSITIVE}),
            "by_kit": suppress({kt: by_kit[kt] for kt in KITS}, {}),
        }
        if total < K_DEFAULT:
            tables = {name: dict.fromkeys(cells) for name, cells in tables.items()}
        return {
            "schema": SCHEMA,
            "k": {"default": K_DEFAULT, "sensitive": K_SENSITIVE},
            "total": total if total >= K_DEFAULT else None,
            **tables,
        }


_repo_instance: DonationRepository | None = None


def get_donation_repo() -> DonationRepository:
    global _repo_instance
    if _repo_instance is None:
        _repo_instance = DonationRepository()
    return _repo_instance
