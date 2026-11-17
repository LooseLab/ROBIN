"""Per-user approval flags (training, report export, etc.)."""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Dict, Mapping, Optional

if TYPE_CHECKING:
    from .store import SecurityStore

TRAINING_RECEIVED_KEY = "training_received"
REPORT_EXPORT_KEY = "report_export"
MINKNOW_REMOTE_CONTROL_KEY = "minknow_remote_control"
ADMIN_USER_APPROVALS_UPDATED_EVENT = "admin.user.approvals_updated"

TRAINING_REQUIRED_MESSAGE = (
    "Access requires training approvals from an administrator on this ROBIN instance."
)
EXPORT_DENIED_MESSAGE = (
    "You need report export approval from an administrator on this ROBIN instance."
)
MINKNOW_REMOTE_CONTROL_DENIED_MESSAGE = (
    "You need MinKNOW remote control approval from an administrator on this ROBIN "
    "instance."
)


@dataclass(frozen=True)
class UserApprovalField:
    key: str
    label: str
    admin_table: bool = True


USER_APPROVAL_FIELDS = (
    UserApprovalField(TRAINING_RECEIVED_KEY, "Training Received"),
    UserApprovalField(REPORT_EXPORT_KEY, "Report Export"),
    UserApprovalField(MINKNOW_REMOTE_CONTROL_KEY, "MinKNOW Remote Control"),
)


def approval_field_labels() -> Dict[str, str]:
    return {field.key: field.label for field in USER_APPROVAL_FIELDS}


def default_approvals() -> Dict[str, bool]:
    return {field.key: False for field in USER_APPROVAL_FIELDS}


def normalize_approvals(
    values: Optional[Mapping[str, Any]],
    *,
    existing: Optional[Mapping[str, bool]] = None,
) -> Dict[str, bool]:
    """Merge and sanitise approval booleans."""
    merged: Dict[str, bool] = dict(existing or default_approvals())
    if not values:
        return merged
    allowed = approval_field_labels()
    for raw_key, raw_value in values.items():
        key = str(raw_key).strip()
        if key not in allowed:
            continue
        merged[key] = bool(raw_value)
    return merged


def parse_approvals_json(raw: Optional[str]) -> Dict[str, bool]:
    if not raw:
        return default_approvals()
    try:
        data = json.loads(raw)
    except Exception:
        return default_approvals()
    if not isinstance(data, dict):
        return default_approvals()
    return normalize_approvals(data)


def approvals_to_json(approvals: Mapping[str, bool]) -> str:
    return json.dumps(dict(normalize_approvals(approvals)), separators=(",", ":"))


def approval_changes(
    previous: Mapping[str, bool],
    updated: Mapping[str, bool],
) -> Dict[str, Dict[str, bool]]:
    """Return only approval keys whose value changed."""
    prev = normalize_approvals(previous)
    new = normalize_approvals(updated)
    changes: Dict[str, Dict[str, bool]] = {}
    for key in approval_field_labels():
        if prev.get(key, False) != new.get(key, False):
            changes[key] = {"from": prev.get(key, False), "to": new.get(key, False)}
    return changes


def approval_audit_details(
    previous: Mapping[str, bool],
    updated: Mapping[str, bool],
    *,
    source: str = "",
    extra: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Build audit-event details for an approval update."""
    prev = normalize_approvals(previous)
    new = normalize_approvals(updated)
    details: Dict[str, Any] = {
        "previous": prev,
        "updated": new,
        "changes": approval_changes(prev, new),
    }
    if source:
        details["source"] = source
    if extra:
        details.update(extra)
    return details


def user_has_approval(store: "SecurityStore", user_id: Optional[int], key: str) -> bool:
    """Return effective approval; active admins always have all approvals."""
    if user_id is None:
        return False
    if store.user_has_role(int(user_id), "admin"):
        user = store.get_user_by_id(int(user_id))
        return bool(user and user.is_active)
    user = store.get_user_by_id(int(user_id))
    if user is None or not user.is_active:
        return False
    return bool(user.approvals.get(key, False))


def effective_approvals(store: "SecurityStore", user_id: Optional[int]) -> Dict[str, bool]:
    """Return effective approval map for display (admins show all granted)."""
    labels = approval_field_labels()
    if user_id is None:
        return {key: False for key in labels}
    if store.user_has_role(int(user_id), "admin"):
        user = store.get_user_by_id(int(user_id))
        granted = bool(user and user.is_active)
        return {key: granted for key in labels}
    user = store.get_user_by_id(int(user_id))
    if user is None or not user.is_active:
        return {key: False for key in labels}
    return {key: bool(user.approvals.get(key, False)) for key in labels}
