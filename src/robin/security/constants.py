from __future__ import annotations

import os
from pathlib import Path

DEFAULT_CONSENT_VERSION = "v1"
CONSENT_VERSION_ENV = "ROBIN_CONSENT_VERSION"


def get_consent_version() -> str:
    """Return the active research-use consent version (env override supported)."""
    raw = os.environ.get(CONSENT_VERSION_ENV, DEFAULT_CONSENT_VERSION)
    version = str(raw or DEFAULT_CONSENT_VERSION).strip()
    return version or DEFAULT_CONSENT_VERSION


def get_security_db_path() -> Path:
    """Return the default SQLite path for security/auth/audit data."""
    if os.name == "nt":
        base = Path(os.environ.get("APPDATA", os.path.expanduser("~")))
    else:
        base = Path(os.environ.get("XDG_CONFIG_HOME", os.path.expanduser("~/.config")))
    return base / "robin" / "security.db"
