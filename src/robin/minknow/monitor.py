"""MinKNOW status helpers for CLI, GUI, and one-shot fetches."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

from robin.minknow.client import MinKnowClient, MinKnowConnectionError
from robin.minknow.config import MinKnowSettings
from robin.minknow.models import PositionStatus, SequencerStatus
from robin.minknow.watch import (
    position_has_active_run,
    position_is_watchable,
    preferred_watch_path,
)


@dataclass(frozen=True)
class MinKnowPollResult:
    """Outcome of a single MinKNOW status poll."""

    status: Optional[SequencerStatus]
    error: Optional[str]


def fetch_sequencer_status(settings: MinKnowSettings) -> MinKnowPollResult:
    """Connect to MinKNOW, fetch status once, and disconnect."""
    if not settings.enabled:
        return MinKnowPollResult(None, "MinKNOW monitoring is disabled.")

    try:
        with MinKnowClient(settings.auth) as client:
            return MinKnowPollResult(client.get_sequencer_status(), None)
    except ImportError as exc:
        return MinKnowPollResult(None, str(exc))
    except MinKnowConnectionError as exc:
        return MinKnowPollResult(None, str(exc))
    except Exception as exc:
        return MinKnowPollResult(None, f"Unexpected error: {exc}")


def position_table_rows(status: SequencerStatus) -> list[dict[str, Any]]:
    """Build NiceGUI table rows from sequencer status."""
    rows: list[dict[str, Any]] = []
    for position in status.positions:
        watch_path = preferred_watch_path(position)
        rows.append(
            {
                "position": position.name,
                "state": position.state,
                "protocol_state": position.protocol_state,
                "sample_id": position.sample_id or "—",
                "protocol_run_id": _short_id(position.protocol_run_id),
                "flow_cell_id": position.flow_cell_id or "—",
                "passed_reads": (
                    str(position.passed_reads)
                    if position.passed_reads is not None
                    else "—"
                ),
                "output_path": position.output_path or "—",
                "watch_path": str(watch_path) if watch_path is not None else "—",
                "can_watch": position_is_watchable(position),
                "can_stop": position_has_active_run(position),
                "protocol_run_id_full": position.protocol_run_id or "",
                "error": position.connection_error or "",
            }
        )
    return rows


def format_poll_summary(
    result: MinKnowPollResult,
    *,
    host: str,
    last_updated: Optional[str] = None,
) -> str:
    """One-line summary for compact UI panels."""
    if result.error:
        return f"MinKNOW ({host}): {result.error}"
    if result.status is None:
        return f"MinKNOW ({host}): no data"

    active = sum(
        1 for position in result.status.positions if position_has_active_run(position)
    )
    parts = [
        f"MinKNOW Core {result.status.core_version}",
        f"{len(result.status.positions)} position(s)",
    ]
    if active:
        parts.append(f"{active} active")
    if last_updated:
        parts.append(f"updated {last_updated}")
    return " · ".join(parts)


def _short_id(value: Optional[str]) -> str:
    if not value:
        return "—"
    if len(value) <= 12:
        return value
    return f"{value[:8]}…"
