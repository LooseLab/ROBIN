from __future__ import annotations

import hashlib
import os
import sqlite3
import time
from dataclasses import dataclass
from typing import Dict, Optional, Tuple


def _now() -> float:
    return time.time()


def _safe_stat(path: str) -> Tuple[int, int, int, float]:
    """
    Returns (st_dev, st_ino, st_size, st_mtime).
    Uses zeros if stat fails.
    """
    try:
        st = os.stat(path)
        return (
            int(getattr(st, "st_dev", 0) or 0),
            int(getattr(st, "st_ino", 0) or 0),
            int(getattr(st, "st_size", 0) or 0),
            float(getattr(st, "st_mtime", 0.0) or 0.0),
        )
    except Exception:
        return 0, 0, 0, 0.0


def file_fingerprint_key(path: str) -> str:
    """
    Stable-ish key across runs. Prefers inode identity when available; falls back to path.

    Note: this is not a content hash (too expensive for BAMs). If a file is replaced in-place
    (same path) without changing mtime/size (rare), it could collide.
    """
    p = os.path.realpath(path or "")
    dev, ino, size, mtime = _safe_stat(p)
    raw = f"v1|{dev}|{ino}|{size}|{mtime:.6f}|{p}"
    return hashlib.blake2b(raw.encode("utf-8", "ignore"), digest_size=16).hexdigest()


@dataclass(frozen=True)
class Stage:
    name: str
    bit: int

    @property
    def mask(self) -> int:
        return 1 << self.bit


DEFAULT_STAGES: Dict[str, Stage] = {
    # General
    "seen": Stage("seen", 0),
    # Pipeline stages (done flags)
    "preprocessing": Stage("preprocessing", 1),
    "bed_conversion": Stage("bed_conversion", 2),
    "mgmt": Stage("mgmt", 3),
    "cnv": Stage("cnv", 4),
    "target": Stage("target", 5),
    "fusion": Stage("fusion", 6),
    "itd": Stage("itd", 17),
    "sturgeon": Stage("sturgeon", 7),
    "nanodx": Stage("nanodx", 8),
    "pannanodx": Stage("pannanodx", 9),
    "random_forest": Stage("random_forest", 10),
    "igv_bam": Stage("igv_bam", 11),
    "snp_analysis": Stage("snp_analysis", 12),
    "target_bam_finalize": Stage("target_bam_finalize", 13),
    "marlin": Stage("marlin", 14),
    "lamprey": Stage("lamprey", 15),
    "tucan": Stage("tucan", 16),
}


class SQLiteStateTracker:
    """
    Persistent bitmask-backed state tracking for files/samples.

    Designed to be used from a single coordinating process/actor (e.g. the Ray Coordinator),
    so we keep the API simple and avoid external locking complexity.
    """

    def __init__(self, db_path: str, stages: Optional[Dict[str, Stage]] = None):
        self.db_path = db_path
        self.stages = stages or DEFAULT_STAGES
        os.makedirs(os.path.dirname(db_path) or ".", exist_ok=True)
        self._conn = sqlite3.connect(
            db_path,
            timeout=30.0,
            isolation_level=None,  # autocommit
            check_same_thread=False,
        )
        self._conn.execute("PRAGMA journal_mode=WAL;")
        self._conn.execute("PRAGMA synchronous=NORMAL;")
        self._conn.execute("PRAGMA busy_timeout=30000;")
        self._init_schema()

    def close(self) -> None:
        try:
            self._conn.close()
        except Exception:
            pass

    def _init_schema(self) -> None:
        self._conn.execute("""
            CREATE TABLE IF NOT EXISTS file_state (
              file_key TEXT PRIMARY KEY,
              path TEXT,
              dev INTEGER,
              ino INTEGER,
              size INTEGER,
              mtime REAL,
              sample_id TEXT,
              flags INTEGER NOT NULL DEFAULT 0,
              first_seen REAL,
              last_seen REAL,
              last_update REAL,
              last_error TEXT
            );
            """)
        self._conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_file_state_sample_id ON file_state(sample_id);"
        )

    def _stage_mask(self, stage_name: str) -> int:
        st = self.stages.get(stage_name)
        if not st:
            return 0
        return st.mask

    def _upsert_seen(
        self, path: str, sample_id: Optional[str] = None, err: Optional[str] = None
    ) -> str:
        p = os.path.realpath(path or "")
        dev, ino, size, mtime = _safe_stat(p)
        fk = file_fingerprint_key(p)
        ts = _now()
        seen_mask = self._stage_mask("seen")
        self._conn.execute("BEGIN;")
        try:
            self._conn.execute(
                """
                INSERT INTO file_state(file_key, path, dev, ino, size, mtime, sample_id, flags, first_seen, last_seen, last_update, last_error)
                VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(file_key) DO UPDATE SET
                  path=excluded.path,
                  dev=excluded.dev,
                  ino=excluded.ino,
                  size=excluded.size,
                  mtime=excluded.mtime,
                  sample_id=COALESCE(excluded.sample_id, file_state.sample_id),
                  flags=(file_state.flags | ?),
                  last_seen=?,
                  last_update=?,
                  last_error=COALESCE(?, file_state.last_error);
                """,
                (
                    fk,
                    p,
                    dev,
                    ino,
                    size,
                    mtime,
                    sample_id,
                    seen_mask,
                    ts,
                    ts,
                    ts,
                    err,
                    seen_mask,
                    ts,
                    ts,
                    err,
                ),
            )
            self._conn.execute("COMMIT;")
        except Exception:
            try:
                self._conn.execute("ROLLBACK;")
            except Exception:
                pass
            raise
        return fk

    def get_flags(self, path: str) -> int:
        fk = file_fingerprint_key(path)
        row = self._conn.execute(
            "SELECT flags FROM file_state WHERE file_key=?;", (fk,)
        ).fetchone()
        return int(row[0]) if row else 0

    def is_done(self, path: str, stage_name: str) -> bool:
        mask = self._stage_mask(stage_name)
        if mask == 0:
            return False
        return (self.get_flags(path) & mask) != 0

    def mark_done(
        self,
        path: str,
        stage_name: str,
        sample_id: Optional[str] = None,
    ) -> None:
        fk = self._upsert_seen(path, sample_id=sample_id)
        mask = self._stage_mask(stage_name)
        if mask == 0:
            return
        ts = _now()
        self._conn.execute(
            """
            UPDATE file_state
            SET flags=(flags | ?), last_update=?, last_error=NULL, sample_id=COALESCE(?, sample_id)
            WHERE file_key=?;
            """,
            (mask, ts, sample_id, fk),
        )

    def mark_failed(
        self,
        path: str,
        stage_name: str,
        error: Optional[str] = None,
        sample_id: Optional[str] = None,
    ) -> None:
        # We still set "seen" and update last_error for post-mortem.
        fk = self._upsert_seen(path, sample_id=sample_id, err=(error or "")[:2000])
        ts = _now()
        # Optional: record failure without consuming a stage bit (keeps bitspace for "done" flags).
        self._conn.execute(
            """
            UPDATE file_state
            SET last_update=?, last_error=COALESCE(?, last_error), sample_id=COALESCE(?, sample_id)
            WHERE file_key=?;
            """,
            (ts, (error or "")[:2000], sample_id, fk),
        )

    def should_skip(
        self,
        path: str,
        job_type: str,
        force: bool = False,
    ) -> bool:
        """
        Returns True if this job_type is already marked done for the given path.
        """
        if force:
            return False
        mask = self._stage_mask(job_type)
        if mask == 0:
            return False
        return self.is_done(path, job_type)
