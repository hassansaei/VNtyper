"""Repository layer for storing and aggregating anonymous data donations."""

from __future__ import annotations

import json
import logging
import sqlite3
import uuid
from datetime import datetime, timezone
from typing import Any

from ..config import settings

logger = logging.getLogger(__name__)

CREATE_TABLE_SQL = """
CREATE TABLE IF NOT EXISTS donations (
    id TEXT PRIMARY KEY,
    run_id TEXT NOT NULL,
    tool_version TEXT NOT NULL,
    sample_hash TEXT NOT NULL,
    decision_files_digest TEXT NOT NULL,
    report_integrity_digest TEXT NOT NULL,
    kit TEXT NOT NULL,
    sequencing_platform TEXT NOT NULL,
    phenotype_hpo TEXT NOT NULL,
    positive_call BOOLEAN NOT NULL,
    confirmation_method TEXT,
    depth_counting_policy TEXT NOT NULL,
    mean_coverage REAL,
    flank_mean_depth REAL,
    sex TEXT,
    collection_month TEXT,
    created_at TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_donations_kit ON donations (kit);
CREATE INDEX IF NOT EXISTS idx_donations_platform ON donations (sequencing_platform);
CREATE INDEX IF NOT EXISTS idx_donations_run_id ON donations (run_id);
"""


class DonationRepository:
    def __init__(self, db_url: str | None = None) -> None:
        self.db_url = db_url or settings.POSTGRES_URL
        self._is_postgres = bool(
            self.db_url and (self.db_url.startswith("postgresql://") or self.db_url.startswith("postgres://"))
        )
        self._init_db()

    def _get_connection(self):
        if self._is_postgres:
            import psycopg

            return psycopg.connect(self.db_url)
        else:
            # Default to SQLite (in memory or local file if url specified)
            sqlite_path = ":memory:"
            if self.db_url and not self.db_url.startswith(("postgres://", "postgresql://")):
                sqlite_path = self.db_url
            conn = sqlite3.connect(sqlite_path)
            conn.row_factory = sqlite3.Row
            return conn

    def _init_db(self) -> None:
        try:
            conn = self._get_connection()
            try:
                cur = conn.cursor()
                try:
                    if self._is_postgres:
                        cur.execute(CREATE_TABLE_SQL)
                    else:
                        cur.executescript(CREATE_TABLE_SQL)
                    conn.commit()
                finally:
                    cur.close()
            finally:
                conn.close()
            logger.info("Donation database initialized successfully.")
        except (sqlite3.Error, OSError, RuntimeError, ValueError, ImportError) as e:
            logger.warning("Could not initialize donation database: %s", e)

    def save_donation(self, data: dict[str, Any]) -> str:
        donation_id = str(uuid.uuid4())
        created_at = datetime.now(timezone.utc).isoformat()
        hpo_json = json.dumps(data.get("phenotype_hpo", []))

        params = (
            donation_id,
            data["run_id"],
            data["tool_version"],
            data["sample_hash"],
            data["decision_files_digest"],
            data["report_integrity_digest"],
            data["kit"],
            data["sequencing_platform"],
            hpo_json,
            bool(data["positive_call"]),
            data.get("confirmation_method"),
            data.get("depth_counting_policy", "vntr_flank_mean_depth"),
            data.get("mean_coverage"),
            data.get("flank_mean_depth"),
            data.get("sex"),
            data.get("collection_month"),
            created_at,
        )

        insert_sql = (
            """
        INSERT INTO donations (
            id, run_id, tool_version, sample_hash, decision_files_digest,
            report_integrity_digest, kit, sequencing_platform, phenotype_hpo,
            positive_call, confirmation_method, depth_counting_policy,
            mean_coverage, flank_mean_depth, sex, collection_month, created_at
        ) VALUES (
            %s, %s, %s, %s, %s,
            %s, %s, %s, %s,
            %s, %s, %s,
            %s, %s, %s, %s, %s
        )
        """
            if self._is_postgres
            else """
        INSERT INTO donations (
            id, run_id, tool_version, sample_hash, decision_files_digest,
            report_integrity_digest, kit, sequencing_platform, phenotype_hpo,
            positive_call, confirmation_method, depth_counting_policy,
            mean_coverage, flank_mean_depth, sex, collection_month, created_at
        ) VALUES (
            ?, ?, ?, ?, ?,
            ?, ?, ?, ?,
            ?, ?, ?,
            ?, ?, ?, ?, ?
        )
        """
        )

        conn = self._get_connection()
        try:
            cur = conn.cursor()
            try:
                cur.execute(insert_sql, params)
                conn.commit()
            finally:
                cur.close()
        finally:
            conn.close()

        return donation_id

    def get_aggregates(self, min_cell_size: int = 5) -> dict[str, Any]:
        """Aggregate donation records enforcing cell-size suppression and negative power checks."""
        conn = self._get_connection()
        try:
            cur = conn.cursor()
            try:
                # Total counts
                cur.execute("SELECT COUNT(*), SUM(CASE WHEN positive_call THEN 1 ELSE 0 END) FROM donations")
                row = cur.fetchone()
                total = row[0] if row else 0
                pos_count = row[1] if row and row[1] is not None else 0
                neg_count = total - pos_count

                # By kit
                cur.execute("""
                    SELECT
                        kit,
                        COUNT(*) as cnt,
                        SUM(CASE WHEN positive_call THEN 1 ELSE 0 END) as pos_cnt,
                        SUM(CASE WHEN NOT positive_call THEN 1 ELSE 0 END) as neg_cnt,
                        AVG(mean_coverage) as avg_cov
                    FROM donations
                    GROUP BY kit
                """)
                kit_rows = cur.fetchall()

                # By platform
                cur.execute("SELECT sequencing_platform, COUNT(*) FROM donations GROUP BY sequencing_platform")
                platform_rows = cur.fetchall()

                # HPO terms
                cur.execute("SELECT phenotype_hpo FROM donations")
                hpo_rows = cur.fetchall()
            finally:
                cur.close()
        finally:
            conn.close()

        by_kit: dict[str, Any] = {}
        suppressed_kit_count = 0
        suppressed_kit_pos = 0
        suppressed_kit_neg = 0

        for r in kit_rows:
            kit_name = r[0]
            cnt = r[1]
            pos = r[2] if r[2] is not None else 0
            neg = r[3] if r[3] is not None else 0
            cov = float(r[4]) if r[4] is not None else None

            if cnt < min_cell_size:
                suppressed_kit_count += cnt
                suppressed_kit_pos += pos
                suppressed_kit_neg += neg
            else:
                by_kit[kit_name] = {
                    "total_samples": cnt,
                    "positive_count": pos,
                    "negative_count": neg,
                    "negative_power_sufficient": (neg >= min_cell_size),
                    "mean_coverage": round(cov, 2) if cov is not None else None,
                }

        if suppressed_kit_count > 0:
            by_kit["Other (<5 samples)"] = {
                "total_samples": suppressed_kit_count,
                "positive_count": suppressed_kit_pos,
                "negative_count": suppressed_kit_neg,
                "negative_power_sufficient": (suppressed_kit_neg >= min_cell_size),
                "mean_coverage": None,
            }

        by_platform: dict[str, int] = {}
        suppressed_plat_count = 0
        for r in platform_rows:
            plat_name = r[0]
            cnt = r[1]
            if cnt < min_cell_size:
                suppressed_plat_count += cnt
            else:
                by_platform[plat_name] = cnt

        if suppressed_plat_count > 0:
            by_platform["Other (<5 samples)"] = suppressed_plat_count

        # Count HPO terms
        hpo_counts: dict[str, int] = {}
        for r in hpo_rows:
            raw = r[0]
            if raw:
                try:
                    terms = json.loads(raw)
                    for t in terms:
                        hpo_counts[t] = hpo_counts.get(t, 0) + 1
                except (json.JSONDecodeError, TypeError):
                    pass

        # Minimum cell size filter for HPO
        filtered_hpo = {k: v for k, v in hpo_counts.items() if v >= min_cell_size}

        return {
            "total_donations": total,
            "positive_count": pos_count,
            "negative_count": neg_count,
            "by_kit": by_kit,
            "by_platform": by_platform,
            "top_hpo_terms": filtered_hpo,
            "minimum_cell_size_threshold": min_cell_size,
        }


# Singleton repository instance
_repo_instance: DonationRepository | None = None


def get_donation_repo() -> DonationRepository:
    global _repo_instance
    if _repo_instance is None:
        _repo_instance = DonationRepository()
    return _repo_instance
