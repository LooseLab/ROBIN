"""Optional dependency helpers for MinKNOW integration."""

from __future__ import annotations

import importlib.util
from typing import Any


def minknow_api_available() -> bool:
    """Return whether ``minknow_api`` can be imported in this environment."""
    return importlib.util.find_spec("minknow_api") is not None


def require_minknow_api() -> Any:
    """Import ``minknow_api`` or raise a clear error."""
    try:
        import minknow_api  # type: ignore import-not-found
    except ImportError as exc:
        raise ImportError(
            "MinKNOW integration requires the minknow_api package. "
            "Install with: pip install 'robin[minknow]'"
        ) from exc
    return minknow_api
