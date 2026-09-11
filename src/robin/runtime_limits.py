"""Process-wide limits for native thread pools used by analysis libraries."""

from __future__ import annotations

import os
from typing import Dict, Optional


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

_BAM_IO_THREADS_MAX = 16


def bam_io_threads_for_preset(
    preset: Optional[str], *, cpu_count: Optional[int] = None
) -> int:
    """BGZF threads for pysam merge/sort/index. Caps at available CPUs."""
    n = max(1, int(cpu_count if cpu_count is not None else (os.cpu_count() or 4)))
    name = (preset or "standard").lower().strip() or "standard"
    if name == "p2i":
        desired = 1
    elif name == "high":
        # One Target actor can use spare cores; leave the mapping conservative.
        desired = max(8, n // 2)
        desired = min(desired, _BAM_IO_THREADS_MAX)
    else:
        desired = 4
    return max(1, min(desired, n))


def resolve_bam_io_threads(
    explicit: Optional[int] = None, *, cpu_count: Optional[int] = None
) -> int:
    """Prefer an explicit value, then ROBIN_BAM_IO_THREADS, then the preset."""
    if explicit is not None:
        try:
            return max(1, int(explicit))
        except (TypeError, ValueError):
            pass
    raw = os.environ.get("ROBIN_BAM_IO_THREADS", "").strip()
    if raw:
        try:
            return max(1, int(raw))
        except ValueError:
            pass
    return bam_io_threads_for_preset(
        os.environ.get("ROBIN_PRESET"), cpu_count=cpu_count
    )


def publish_workflow_preset(preset: Optional[str]) -> int:
    """Export preset and default BAM I/O threads for Ray workers."""
    name = (preset or "standard").lower().strip() or "standard"
    os.environ["ROBIN_PRESET"] = name
    threads = bam_io_threads_for_preset(name)
    os.environ.setdefault("ROBIN_BAM_IO_THREADS", str(threads))
    return threads


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
