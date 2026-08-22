"""Process-wide limits for native thread pools used by analysis libraries."""

from __future__ import annotations

import os
from typing import Dict

_NATIVE_THREAD_ENV_VARS = (
    "POLARS_MAX_THREADS",
    "RAYON_NUM_THREADS",
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS",
)


def configure_native_thread_limits() -> Dict[str, str]:
    """Cap native pools before importing Polars, NumPy, or analysis modules."""
    raw_value = os.environ.get("ROBIN_NATIVE_THREADS", "1")
    try:
        thread_count = max(1, int(raw_value))
    except (TypeError, ValueError):
        thread_count = 1

    value = str(thread_count)
    configured: Dict[str, str] = {}
    for name in _NATIVE_THREAD_ENV_VARS:
        os.environ.setdefault(name, value)
        configured[name] = os.environ[name]
    return configured
