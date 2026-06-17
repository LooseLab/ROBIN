from __future__ import annotations

import json
import sqlite3
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from .constants import get_security_db_path
from .models import User, UserPublic


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


class SecurityStore:
    """SQLite-backed store for users, roles, consent, and audit events."""

    def __init__(self, db_path: Optional[Path] = None):
        self.db_path = Path(db_path) if db_path else get_security_db_path()
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._lock = threading.RLock()
        self._conn = sqlite3.connect(
            str(self.db_path),
            timeout=30.0,
            isolation_level=None,
            check_same_thread=False,
        )
        self._conn.row_factory = sqlite3.Row
        self._conn.execute("PRAGMA journal_mode=WAL;")
        self._conn.execute("PRAGMA synchronous=NORMAL;")
        self._conn.execute("PRAGMA busy_timeout=30000;")
        self._init_schema()
        self._ensure_default_roles()

    def _init_schema(self) -> None:
        with self._lock:
            self._conn.executescript(
                """
                CREATE TABLE IF NOT EXISTS users (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    username TEXT NOT NULL UNIQUE,
                    password_hash TEXT NOT NULL,
                    is_active INTEGER NOT NULL DEFAULT 1,
                    created_at TEXT NOT NULL,
                    last_login_at TEXT
                );

                CREATE TABLE IF NOT EXISTS roles (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    name TEXT NOT NULL UNIQUE
                );

                CREATE TABLE IF NOT EXISTS user_roles (
                    user_id INTEGER NOT NULL,
                    role_id INTEGER NOT NULL,
                    PRIMARY KEY (user_id, role_id),
                    FOREIGN KEY(user_id) REFERENCES users(id),
                    FOREIGN KEY(role_id) REFERENCES roles(id)
                );

                CREATE TABLE IF NOT EXISTS consents (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    user_id INTEGER NOT NULL,
                    consent_version TEXT NOT NULL,
                    agreed_at TEXT NOT NULL,
                    ip TEXT,
                    user_agent TEXT,
                    session_id TEXT,
                    FOREIGN KEY(user_id) REFERENCES users(id)
                );

                CREATE TABLE IF NOT EXISTS audit_events (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    occurred_at TEXT NOT NULL,
                    user_id INTEGER,
                    event_type TEXT NOT NULL,
                    target_type TEXT,
                    target_id TEXT,
                    details_json TEXT,
                    ip TEXT,
                    user_agent TEXT,
                    session_id TEXT,
                    request_id TEXT,
                    result TEXT NOT NULL,
                    error_code TEXT,
                    FOREIGN KEY(user_id) REFERENCES users(id)
                );

                CREATE INDEX IF NOT EXISTS idx_audit_events_occurred_at ON audit_events(occurred_at);
                CREATE INDEX IF NOT EXISTS idx_audit_events_user_id ON audit_events(user_id);
                CREATE INDEX IF NOT EXISTS idx_audit_events_event_type ON audit_events(event_type);
                CREATE INDEX IF NOT EXISTS idx_audit_events_target_id ON audit_events(target_id);
                CREATE INDEX IF NOT EXISTS idx_consents_user_version ON consents(user_id, consent_version);
                """
            )

    def _ensure_default_roles(self) -> None:
        with self._lock:
            self._conn.execute(
                "INSERT OR IGNORE INTO roles(name) VALUES (?)",
                ("admin",),
            )
            self._conn.execute(
                "INSERT OR IGNORE INTO roles(name) VALUES (?)",
                ("user",),
            )

    def has_users(self) -> bool:
        with self._lock:
            row = self._conn.execute("SELECT COUNT(*) AS c FROM users").fetchone()
            return bool(row and int(row["c"]) > 0)

    def create_user(self, username: str, password_hash: str, *, is_active: bool = True) -> int:
        created_at = utc_now_iso()
        with self._lock:
            cur = self._conn.execute(
                "INSERT INTO users(username, password_hash, is_active, created_at) VALUES (?, ?, ?, ?)",
                (username.strip(), password_hash, 1 if is_active else 0, created_at),
            )
            return int(cur.lastrowid)

    def assign_role(self, user_id: int, role_name: str) -> None:
        with self._lock:
            role = self._conn.execute(
                "SELECT id FROM roles WHERE name = ?",
                (role_name,),
            ).fetchone()
            if role is None:
                raise ValueError(f"Unknown role: {role_name}")
            self._conn.execute(
                "INSERT OR IGNORE INTO user_roles(user_id, role_id) VALUES (?, ?)",
                (user_id, int(role["id"])),
            )

    def get_user_by_username(self, username: str) -> Optional[User]:
        with self._lock:
            row = self._conn.execute(
                "SELECT id, username, password_hash, is_active, created_at, last_login_at FROM users WHERE username = ?",
                (username.strip(),),
            ).fetchone()
            if row is None:
                return None
            return User(
                id=int(row["id"]),
                username=str(row["username"]),
                password_hash=str(row["password_hash"]),
                is_active=bool(row["is_active"]),
                created_at=str(row["created_at"]),
                last_login_at=str(row["last_login_at"]) if row["last_login_at"] else None,
            )

    def get_user_by_id(self, user_id: int) -> Optional[User]:
        with self._lock:
            row = self._conn.execute(
                "SELECT id, username, password_hash, is_active, created_at, last_login_at FROM users WHERE id = ?",
                (int(user_id),),
            ).fetchone()
            if row is None:
                return None
            return User(
                id=int(row["id"]),
                username=str(row["username"]),
                password_hash=str(row["password_hash"]),
                is_active=bool(row["is_active"]),
                created_at=str(row["created_at"]),
                last_login_at=str(row["last_login_at"]) if row["last_login_at"] else None,
            )

    def get_user_public(self, user_id: int) -> Optional[UserPublic]:
        with self._lock:
            row = self._conn.execute(
                "SELECT id, username, is_active, created_at, last_login_at FROM users WHERE id = ?",
                (int(user_id),),
            ).fetchone()
            if row is None:
                return None
            return UserPublic(
                id=int(row["id"]),
                username=str(row["username"]),
                is_active=bool(row["is_active"]),
                created_at=str(row["created_at"]),
                last_login_at=str(row["last_login_at"]) if row["last_login_at"] else None,
            )

    def set_user_password_hash(self, username: str, password_hash: str) -> bool:
        with self._lock:
            cur = self._conn.execute(
                "UPDATE users SET password_hash = ? WHERE username = ?",
                (password_hash, username.strip()),
            )
            return int(cur.rowcount) > 0

    def set_user_active(self, username: str, is_active: bool) -> bool:
        with self._lock:
            cur = self._conn.execute(
                "UPDATE users SET is_active = ? WHERE username = ?",
                (1 if is_active else 0, username.strip()),
            )
            return int(cur.rowcount) > 0

    def list_users(self) -> List[UserPublic]:
        with self._lock:
            rows = self._conn.execute(
                "SELECT id, username, is_active, created_at, last_login_at FROM users ORDER BY username"
            ).fetchall()
            return [
                UserPublic(
                    id=int(row["id"]),
                    username=str(row["username"]),
                    is_active=bool(row["is_active"]),
                    created_at=str(row["created_at"]),
                    last_login_at=str(row["last_login_at"]) if row["last_login_at"] else None,
                )
                for row in rows
            ]

    def user_has_role(self, user_id: int, role_name: str) -> bool:
        with self._lock:
            row = self._conn.execute(
                """
                SELECT 1
                FROM user_roles ur
                INNER JOIN roles r ON r.id = ur.role_id
                WHERE ur.user_id = ? AND r.name = ?
                LIMIT 1
                """,
                (int(user_id), role_name),
            ).fetchone()
            return row is not None

    def revoke_role(self, user_id: int, role_name: str) -> bool:
        with self._lock:
            role = self._conn.execute(
                "SELECT id FROM roles WHERE name = ?",
                (role_name,),
            ).fetchone()
            if role is None:
                return False
            cur = self._conn.execute(
                "DELETE FROM user_roles WHERE user_id = ? AND role_id = ?",
                (int(user_id), int(role["id"])),
            )
            return int(cur.rowcount) > 0

    def count_active_admins(self) -> int:
        with self._lock:
            row = self._conn.execute(
                """
                SELECT COUNT(*) AS c
                FROM users u
                INNER JOIN user_roles ur ON ur.user_id = u.id
                INNER JOIN roles r ON r.id = ur.role_id
                WHERE u.is_active = 1 AND r.name = 'admin'
                """
            ).fetchone()
            return int(row["c"]) if row else 0

    def set_last_login(self, user_id: int) -> None:
        with self._lock:
            self._conn.execute(
                "UPDATE users SET last_login_at = ? WHERE id = ?",
                (utc_now_iso(), int(user_id)),
            )

    def get_user_roles(self, user_id: int) -> List[str]:
        with self._lock:
            rows = self._conn.execute(
                """
                SELECT r.name
                FROM roles r
                INNER JOIN user_roles ur ON ur.role_id = r.id
                WHERE ur.user_id = ?
                ORDER BY r.name
                """,
                (int(user_id),),
            ).fetchall()
            return [str(r["name"]) for r in rows]

    def has_consent(self, user_id: int, consent_version: str) -> bool:
        with self._lock:
            row = self._conn.execute(
                """
                SELECT 1
                FROM consents
                WHERE user_id = ? AND consent_version = ?
                LIMIT 1
                """,
                (int(user_id), consent_version),
            ).fetchone()
            return row is not None

    def any_active_admin_has_consent(self, consent_version: str) -> bool:
        """Return True when an active admin has accepted the consent version."""
        for user in self.list_users():
            if not user.is_active:
                continue
            if self.user_has_role(user.id, "admin") and self.has_consent(
                user.id, consent_version
            ):
                return True
        return False

    def record_consent(
        self,
        user_id: int,
        consent_version: str,
        *,
        ip: str = "",
        user_agent: str = "",
        session_id: str = "",
    ) -> None:
        with self._lock:
            self._conn.execute(
                """
                INSERT INTO consents(user_id, consent_version, agreed_at, ip, user_agent, session_id)
                VALUES (?, ?, ?, ?, ?, ?)
                """,
                (int(user_id), consent_version, utc_now_iso(), ip, user_agent, session_id),
            )

    def list_consent_status(self, consent_version: str) -> List[Dict[str, Any]]:
        """Return per-user consent status for the given consent version."""
        with self._lock:
            rows = self._conn.execute(
                """
                SELECT
                    u.id AS user_id,
                    u.username,
                    u.is_active,
                    (
                        SELECT c.agreed_at
                        FROM consents c
                        WHERE c.user_id = u.id AND c.consent_version = ?
                        ORDER BY c.id DESC
                        LIMIT 1
                    ) AS agreed_at
                FROM users u
                ORDER BY u.username
                """,
                (consent_version,),
            ).fetchall()
        out: List[Dict[str, Any]] = []
        for row in rows:
            agreed_at = row["agreed_at"]
            out.append(
                {
                    "user_id": int(row["user_id"]),
                    "username": str(row["username"]),
                    "is_active": bool(row["is_active"]),
                    "consent_version": consent_version,
                    "has_consent": agreed_at is not None,
                    "agreed_at": str(agreed_at) if agreed_at else None,
                }
            )
        return out

    def append_audit_event(
        self,
        *,
        event_type: str,
        result: str,
        user_id: Optional[int] = None,
        target_type: str = "",
        target_id: str = "",
        details: Optional[Dict[str, Any]] = None,
        ip: str = "",
        user_agent: str = "",
        session_id: str = "",
        request_id: str = "",
        error_code: str = "",
    ) -> None:
        details_json = json.dumps(details or {}, separators=(",", ":"), ensure_ascii=True)
        with self._lock:
            self._conn.execute(
                """
                INSERT INTO audit_events(
                    occurred_at, user_id, event_type, target_type, target_id, details_json,
                    ip, user_agent, session_id, request_id, result, error_code
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    utc_now_iso(),
                    user_id,
                    event_type,
                    target_type or None,
                    target_id or None,
                    details_json,
                    ip or None,
                    user_agent or None,
                    session_id or None,
                    request_id or None,
                    result,
                    error_code or None,
                ),
            )

    def query_audit_events(
        self,
        *,
        username: str = "",
        event_type: str = "",
        target_type: str = "",
        target_id: str = "",
        sample_id: str = "",
        from_ts: str = "",
        to_ts: str = "",
        limit: int = 200,
    ) -> List[Dict[str, Any]]:
        clauses = []
        params: List[Any] = []
        if username:
            clauses.append("u.username = ?")
            params.append(username.strip())
        if event_type:
            clauses.append("a.event_type = ?")
            params.append(event_type.strip())
        if target_type:
            clauses.append("a.target_type = ?")
            params.append(target_type.strip())
        if target_id:
            clauses.append("a.target_id = ?")
            params.append(target_id.strip())
        if sample_id:
            sid = sample_id.strip()
            sample_json = json.dumps(sid)
            clauses.append(
                "((a.target_type = 'sample' AND a.target_id = ?) "
                "OR instr(a.details_json, ?) > 0)"
            )
            params.extend([sid, sample_json])
        if from_ts:
            clauses.append("a.occurred_at >= ?")
            params.append(from_ts.strip())
        if to_ts:
            clauses.append("a.occurred_at <= ?")
            params.append(to_ts.strip())
        where_sql = f"WHERE {' AND '.join(clauses)}" if clauses else ""
        sql = f"""
            SELECT
                a.id,
                a.occurred_at,
                a.user_id,
                COALESCE(u.username, '') AS username,
                a.event_type,
                a.target_type,
                a.target_id,
                a.details_json,
                a.ip,
                a.user_agent,
                a.session_id,
                a.request_id,
                a.result,
                a.error_code
            FROM audit_events a
            LEFT JOIN users u ON u.id = a.user_id
            {where_sql}
            ORDER BY a.id DESC
            LIMIT ?
        """
        params.append(max(1, int(limit)))
        with self._lock:
            rows = self._conn.execute(sql, tuple(params)).fetchall()
        out: List[Dict[str, Any]] = []
        for row in rows:
            details_raw = str(row["details_json"] or "{}")
            try:
                details = json.loads(details_raw)
            except Exception:
                details = {"raw": details_raw}
            out.append(
                {
                    "id": int(row["id"]),
                    "occurred_at": str(row["occurred_at"]),
                    "user_id": int(row["user_id"]) if row["user_id"] is not None else None,
                    "username": str(row["username"] or ""),
                    "event_type": str(row["event_type"] or ""),
                    "target_type": str(row["target_type"] or ""),
                    "target_id": str(row["target_id"] or ""),
                    "details": details,
                    "ip": str(row["ip"] or ""),
                    "user_agent": str(row["user_agent"] or ""),
                    "session_id": str(row["session_id"] or ""),
                    "request_id": str(row["request_id"] or ""),
                    "result": str(row["result"] or ""),
                    "error_code": str(row["error_code"] or ""),
                }
            )
        return out
