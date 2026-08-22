"""Rules for when MNP-Flex may start without explicit per-sample user action."""

from __future__ import annotations

import csv
import os
import time
from pathlib import Path
from typing import Any, Mapping, Optional

# Matches gui_launcher.GUILauncher.completion_timeout_seconds (Live → Complete).
DEFAULT_MNPFLEX_IDLE_SECONDS = 15 * 60


def mnpflex_idle_seconds() -> int:
    """Seconds with no new sample data before automatic MNP-Flex may run."""
    raw = os.getenv("MNPFLEX_IDLE_SECONDS", str(DEFAULT_MNPFLEX_IDLE_SECONDS))
    try:
        return max(0, int(raw))
    except ValueError:
        return DEFAULT_MNPFLEX_IDLE_SECONDS


def _int_field(data: Mapping[str, Any], *keys: str, default: int = 0) -> int:
    for key in keys:
        if key in data and data[key] not in (None, ""):
            try:
                return int(data[key])
            except (TypeError, ValueError):
                continue
    return default


def _float_field(data: Mapping[str, Any], *keys: str, default: float = 0.0) -> float:
    for key in keys:
        if key in data and data[key] not in (None, ""):
            try:
                return float(data[key])
            except (TypeError, ValueError):
                continue
    return default


def sample_workflow_jobs_complete(overview: Mapping[str, Any]) -> bool:
    active = _int_field(overview, "active_jobs", "samples_overview_active_jobs")
    pending = _int_field(overview, "pending_jobs", "samples_overview_pending_jobs")
    total = _int_field(overview, "total_jobs", "samples_overview_total_jobs")
    completed = _int_field(
        overview, "completed_jobs", "samples_overview_completed_jobs"
    )
    return total > 0 and completed >= total and active == 0 and pending == 0


def sample_data_last_seen(overview: Mapping[str, Any]) -> float:
    return _float_field(
        overview,
        "_last_seen_raw",
        "samples_overview_last_seen",
        "last_seen",
    )


def sample_data_idle_long_enough(
    overview: Mapping[str, Any],
    *,
    idle_seconds: Optional[int] = None,
    now: Optional[float] = None,
) -> bool:
    idle = mnpflex_idle_seconds() if idle_seconds is None else idle_seconds
    if idle <= 0:
        return True
    last_seen = sample_data_last_seen(overview)
    if last_seen <= 0:
        return False
    current = time.time() if now is None else now
    return (current - last_seen) >= idle


def sample_ready_for_mnpflex_auto_run(
    overview: Mapping[str, Any],
    *,
    idle_seconds: Optional[int] = None,
    now: Optional[float] = None,
) -> bool:
    return sample_workflow_jobs_complete(overview) and sample_data_idle_long_enough(
        overview,
        idle_seconds=idle_seconds,
        now=now,
    )


def read_master_csv_overview_row(sample_dir: Path) -> Optional[dict[str, Any]]:
    master_csv = sample_dir / "master.csv"
    if not master_csv.exists():
        return None
    try:
        with master_csv.open("r", newline="") as fh:
            reader = csv.DictReader(fh)
            return next(reader, None)
    except Exception:
        return None


def sample_ready_for_mnpflex_auto_run_from_dir(
    sample_dir: Path,
    *,
    idle_seconds: Optional[int] = None,
    now: Optional[float] = None,
) -> bool:
    row = read_master_csv_overview_row(sample_dir)
    if not row:
        return False
    return sample_ready_for_mnpflex_auto_run(row, idle_seconds=idle_seconds, now=now)
