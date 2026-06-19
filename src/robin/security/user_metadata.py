"""User profile metadata keys and normalisation helpers."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from typing import Any, Dict, Mapping, Optional

EMAIL_KEY = "email"
CLINICAL_ROLE_KEY = "clinical_role"
NOTES_KEY = "notes"

_EMAIL_RE = re.compile(r"^[^@\s]+@[^@\s]+\.[^@\s]+$")


@dataclass(frozen=True)
class UserMetadataField:
    key: str
    label: str
    admin_table: bool = True


USER_METADATA_FIELDS = (
    UserMetadataField(EMAIL_KEY, "Email"),
    UserMetadataField(CLINICAL_ROLE_KEY, "Clinical role"),
    UserMetadataField(NOTES_KEY, "Notes", admin_table=False),
)


def metadata_field_labels() -> Dict[str, str]:
    return {field.key: field.label for field in USER_METADATA_FIELDS}


def normalize_metadata(
    values: Optional[Mapping[str, Any]],
    *,
    existing: Optional[Mapping[str, str]] = None,
) -> Dict[str, str]:
    """Merge and sanitise metadata values (string keys/values, omit blanks)."""
    merged: Dict[str, str] = dict(existing or {})
    if not values:
        return merged
    allowed = metadata_field_labels()
    for raw_key, raw_value in values.items():
        key = str(raw_key).strip()
        if key not in allowed:
            continue
        if raw_value is None:
            merged.pop(key, None)
            continue
        value = str(raw_value).strip()
        if not value:
            merged.pop(key, None)
            continue
        if key == EMAIL_KEY and not _EMAIL_RE.match(value):
            raise ValueError(f"Invalid email address: {value}")
        merged[key] = value
    return merged


def parse_metadata_json(raw: Optional[str]) -> Dict[str, str]:
    if not raw:
        return {}
    try:
        data = json.loads(raw)
    except Exception:
        return {}
    if not isinstance(data, dict):
        return {}
    return normalize_metadata(data)


def metadata_to_json(metadata: Mapping[str, str]) -> str:
    return json.dumps(dict(metadata), separators=(",", ":"), ensure_ascii=True)
