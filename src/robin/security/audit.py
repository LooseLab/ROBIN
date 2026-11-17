from __future__ import annotations

from typing import Any, Dict, Optional

from .store import SecurityStore


class AuditService:
    def __init__(self, store: SecurityStore):
        self.store = store

    def log_event(
        self,
        *,
        event_type: str,
        result: str = "success",
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
        self.store.append_audit_event(
            event_type=event_type,
            result=result,
            user_id=user_id,
            target_type=target_type,
            target_id=target_id,
            details=details,
            ip=ip,
            user_agent=user_agent,
            session_id=session_id,
            request_id=request_id,
            error_code=error_code,
        )
