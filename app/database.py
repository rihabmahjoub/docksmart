"""
Lightweight job persistence.

We use plain sqlite3 (stdlib) rather than an ORM: DockSmart's job table
is a single flat structure, and this keeps the deployment footprint
minimal for free-tier hosts (PythonAnywhere / Render) where every extra
dependency is friction.
"""
from __future__ import annotations

import json
import sqlite3
from contextlib import contextmanager
from datetime import datetime, timezone
from typing import Any, Iterator, Optional

from app.config import settings

SCHEMA = """
CREATE TABLE IF NOT EXISTS jobs (
    job_id          TEXT PRIMARY KEY,
    status          TEXT NOT NULL DEFAULT 'queued',
    stage           TEXT NOT NULL DEFAULT 'created',
    created_at      TEXT NOT NULL,
    updated_at      TEXT NOT NULL,
    receptor_source TEXT,
    ligand_source   TEXT,
    params_json     TEXT,
    result_json     TEXT,
    error_message   TEXT
);
"""


@contextmanager
def get_conn() -> Iterator[sqlite3.Connection]:
    conn = sqlite3.connect(settings.DB_PATH, timeout=30, check_same_thread=False)
    conn.row_factory = sqlite3.Row
    try:
        yield conn
        conn.commit()
    finally:
        conn.close()


def init_db() -> None:
    with get_conn() as conn:
        conn.execute(SCHEMA)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def create_job(job_id: str, receptor_source: str, ligand_source: str, params: dict) -> None:
    with get_conn() as conn:
        conn.execute(
            """INSERT INTO jobs (job_id, status, stage, created_at, updated_at,
                                  receptor_source, ligand_source, params_json)
               VALUES (?, 'queued', 'created', ?, ?, ?, ?, ?)""",
            (job_id, _now(), _now(), receptor_source, ligand_source, json.dumps(params)),
        )


def update_job(
    job_id: str,
    *,
    status: Optional[str] = None,
    stage: Optional[str] = None,
    result: Optional[dict] = None,
    error: Optional[str] = None,
) -> None:
    fields, values = ["updated_at = ?"], [_now()]
    if status is not None:
        fields.append("status = ?")
        values.append(status)
    if stage is not None:
        fields.append("stage = ?")
        values.append(stage)
    if result is not None:
        fields.append("result_json = ?")
        values.append(json.dumps(result))
    if error is not None:
        fields.append("error_message = ?")
        values.append(error)
    values.append(job_id)
    with get_conn() as conn:
        conn.execute(f"UPDATE jobs SET {', '.join(fields)} WHERE job_id = ?", values)


def get_job(job_id: str) -> Optional[dict[str, Any]]:
    with get_conn() as conn:
        row = conn.execute("SELECT * FROM jobs WHERE job_id = ?", (job_id,)).fetchone()
    if row is None:
        return None
    job = dict(row)
    job["params"] = json.loads(job.pop("params_json") or "{}")
    job["result"] = json.loads(job.pop("result_json") or "null")
    return job


def list_recent_jobs(limit: int = 20) -> list[dict[str, Any]]:
    with get_conn() as conn:
        rows = conn.execute(
            "SELECT job_id, status, stage, created_at FROM jobs ORDER BY created_at DESC LIMIT ?",
            (limit,),
        ).fetchall()
    return [dict(r) for r in rows]
