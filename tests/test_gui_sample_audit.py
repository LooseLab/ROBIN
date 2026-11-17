"""Tests for per-sample audit dialog helpers."""

from __future__ import annotations

from pathlib import Path

from robin.gui.components.sample_audit import _audit_rows_for_sample, _export_sample_audit_csv
from robin.security.store import SecurityStore


class _LauncherStub:
    def __init__(self, store: SecurityStore) -> None:
        self.security_store = store
        self.audit_calls: list[dict] = []

    def _get_current_user_id(self):
        return 1

    def _audit_log(self, **kwargs):
        self.audit_calls.append(kwargs)


def _store(tmp_path: Path) -> SecurityStore:
    return SecurityStore(db_path=tmp_path / "security.db")


def test_sample_audit_rows_and_export(tmp_path: Path) -> None:
    store = _store(tmp_path)
    user_id = store.create_user("alice", "hash")
    store.append_audit_event(
        event_type="sample.viewed",
        result="success",
        user_id=user_id,
        target_type="sample",
        target_id="Sample_103",
    )
    store.append_audit_event(
        event_type="run.started",
        result="success",
        user_id=user_id,
        target_type="sample",
        target_id="Sample_104",
        details={"run_type": "mnpflex"},
    )

    launcher = _LauncherStub(store)
    rows = _audit_rows_for_sample(launcher, "Sample_103", limit=10)
    assert len(rows) == 1
    assert rows[0]["event_type"] == "sample.viewed"

    payload = _export_sample_audit_csv(launcher, "Sample_103")
    assert b"sample.viewed" in payload
    assert launcher.audit_calls
    assert launcher.audit_calls[0]["event_type"] == "sample.audit.exported"
