"""Sample ID generation helpers for MinKNOW runs."""

from __future__ import annotations

import hashlib


def generate_sample_id_md5(
    test_id: str,
    *,
    first_name: str = "",
    last_name: str = "",
    date_of_birth: str = "",
) -> str:
    """Return the ROBIN MD5 sample ID used by the Sample ID generator page."""
    parts = [
        test_id.strip(),
        first_name.strip(),
        last_name.strip(),
        date_of_birth.strip(),
    ]
    if not parts[0]:
        raise ValueError("Test ID is required")
    payload = "|".join(parts)
    return hashlib.md5(payload.encode("utf-8")).hexdigest()
