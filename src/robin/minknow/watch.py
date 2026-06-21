"""Watch MinKNOW run output directories in the ROBIN workflow."""

from __future__ import annotations

import logging
import threading
from pathlib import Path
from typing import Optional

from robin.minknow.models import PositionStatus, SequencerStatus

LOGGER = logging.getLogger(__name__)

_INACTIVE_PROTOCOL_STATES = frozenset(
    {
        "no_protocol_state",
        "protocol_finished",
        "protocol_finished_successfully",
        "protocol_finished_failed",
        "protocol_completed",
        "protocol_stopped",
    }
)


def position_has_active_run(position: PositionStatus) -> bool:
    """Return whether a sequencing protocol run is in progress on this position."""
    if not position.sample_id or not position.protocol_run_id:
        return False
    if position.protocol_state in _INACTIVE_PROTOCOL_STATES:
        return False
    return True


def position_is_watchable(position: PositionStatus) -> bool:
    """Return whether a position has an active run worth adding to the workflow watch list."""
    if not position_has_active_run(position):
        return False
    return bool(position.output_path or position.output_reads_path)


def resolve_watch_path(position: PositionStatus) -> tuple[Optional[Path], str]:
    """Pick the best directory to pass to ``add_watch_path``."""
    if not position.sample_id:
        return None, "No sample ID on this run."

    candidates = _candidate_paths(position)
    if not candidates:
        return None, "No output path reported for this run."

    ranked = sorted(candidates, key=lambda path: _candidate_rank(path, position))
    for path in ranked:
        if path.is_dir():
            return path.resolve(), ""

    best = ranked[0]
    if not best.exists():
        return None, f"Output path does not exist yet: {best}"
    if not best.is_dir():
        return None, f"Output path is not a directory: {best}"
    return best.resolve(), ""


def preferred_watch_path(position: PositionStatus) -> Optional[Path]:
    """Return the path that would be watched for an active run."""
    if not position_has_active_run(position):
        return None
    path, _hint = resolve_watch_path(position)
    if path is not None:
        return path
    candidates = _candidate_paths(position)
    return candidates[0] if candidates else None


def watch_position_run(position: PositionStatus) -> tuple[bool, str]:
    """Add the resolved output directory for ``position`` to the workflow watch list."""
    path, hint = resolve_watch_path(position)
    if path is None:
        return False, hint

    try:
        from robin.workflow_ray import add_watch_path
    except ImportError as exc:
        return False, f"Workflow watch integration unavailable: {exc}"

    return add_watch_path(str(path))


def watch_active_runs(status: SequencerStatus) -> list[tuple[str, bool, str]]:
    """Watch every watchable position in ``status`` once."""
    results: list[tuple[str, bool, str]] = []
    for position in status.positions:
        if not position_is_watchable(position):
            continue
        success, message = watch_position_run(position)
        results.append((position.name, success, message))
    return results


class AutoWatchTracker:
    """Track protocol runs that have already been submitted to ``add_watch_path``."""

    def __init__(self) -> None:
        self._watched_runs: set[str] = set()
        self._lock = threading.Lock()

    def reset(self) -> None:
        with self._lock:
            self._watched_runs.clear()

    def process(self, status: SequencerStatus) -> list[tuple[str, bool, str]]:
        """Watch newly detected runs; returns actions taken this update."""
        actions: list[tuple[str, bool, str]] = []
        for position in status.positions:
            if not position_is_watchable(position):
                continue
            run_id = (position.protocol_run_id or "").strip()
            if not run_id:
                continue
            with self._lock:
                if run_id in self._watched_runs:
                    continue

            success, message = watch_position_run(position)
            if success or _should_mark_watched(message):
                with self._lock:
                    self._watched_runs.add(run_id)
            if success:
                actions.append((position.name, True, message))
            elif not _is_retryable_watch_failure(message):
                LOGGER.warning(
                    "Auto-watch failed for %s (%s): %s",
                    position.name,
                    run_id,
                    message,
                )
        return actions


_trackers_lock = threading.Lock()
_trackers: dict[str, AutoWatchTracker] = {}


def get_auto_watch_tracker(host_key: str) -> AutoWatchTracker:
    """Return a shared auto-watch tracker for a MinKNOW host key."""
    with _trackers_lock:
        tracker = _trackers.get(host_key)
        if tracker is None:
            tracker = AutoWatchTracker()
            _trackers[host_key] = tracker
        return tracker


def process_auto_watch(status: SequencerStatus) -> list[tuple[str, bool, str]]:
    """Run auto-watch for a sequencer status update."""
    host_key = f"{status.host}:{status.port}"
    return get_auto_watch_tracker(host_key).process(status)


def _candidate_paths(position: PositionStatus) -> list[Path]:
    candidates: list[Path] = []
    seen: set[str] = set()
    for raw in (position.output_path, position.output_reads_path):
        if not raw or not str(raw).strip() or raw == "—":
            continue
        path = Path(str(raw).strip()).expanduser()
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        candidates.append(path)
    return candidates


def _candidate_rank(path: Path, position: PositionStatus) -> tuple[int, int, int]:
    sample_id = position.sample_id or ""
    sample_specific = int(sample_id and sample_id in str(path))
    exists = int(path.is_dir())
    has_bams = int(_directory_has_bams(path)) if exists else 0
    return (-has_bams, -sample_specific, -exists)


def _directory_has_bams(path: Path) -> bool:
    try:
        if any(path.glob("*.bam")):
            return True
        for child in path.iterdir():
            if child.is_dir() and any(child.glob("*.bam")):
                return True
    except OSError:
        return False
    return False


def _should_mark_watched(message: str) -> bool:
    lowered = message.lower()
    return "already watched" in lowered


def _is_retryable_watch_failure(message: str) -> bool:
    lowered = message.lower()
    return "does not exist" in lowered or "not exist yet" in lowered
