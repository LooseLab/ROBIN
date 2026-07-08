"""Shared helpers for mapping MinKNOW protobuf messages to ROBIN models."""

from __future__ import annotations

from dataclasses import replace
from typing import Any, Optional

from robin.minknow.models import PositionStatus

# Terminal manager protocol states — clear cached run metadata when seen.
_FINISHED_PROTOCOL_STATES = frozenset(
    {
        "protocol_finished",
        "protocol_finished_successfully",
        "protocol_finished_failed",
        "protocol_completed",
        "protocol_stopped",
    }
)

_TERMINAL_ACQUISITION_STATES = frozenset(
    {
        "acquisition_completed",
    }
)

_ACTIVE_ACQUISITION_STATES = frozenset(
    {
        "acquisition_starting",
        "acquisition_running",
        "acquisition_finishing",
    }
)


def position_status_from_description(
    description: Any,
    *,
    host: str,
    credentials: Any,
) -> PositionStatus:
    """Build a ``PositionStatus`` from a manager ``FlowCellPosition`` message."""
    from minknow_api.manager import FlowCellPosition

    position = FlowCellPosition(description, host, credentials)
    return PositionStatus(
        name=position.name,
        state=position.state,
        protocol_state=position.protocol_state,
        running=bool(position.running),
        device_type=getattr(position, "device_type", None),
    )


def merge_position_description(
    current: Optional[PositionStatus],
    incoming: PositionStatus,
) -> PositionStatus:
    """Merge a manager position update into an existing snapshot.

    Manager ``watch_flow_cell_positions`` events only carry high-level fields.
    Activity-stream metadata (sample ID, run ID, output paths, yields) must be
    preserved across manager ``changes`` updates until the run actually ends.
    """
    if current is None or current.name != incoming.name:
        return incoming

    merged = replace(
        current,
        state=incoming.state,
        protocol_state=incoming.protocol_state,
        running=incoming.running,
        device_type=incoming.device_type,
        connection_error=None,
    )

    # Manager ``running`` can flicker during mux/protocol transitions while an
    # acquisition is still active. Only clear cached run metadata when the
    # manager reports a terminal protocol state (or acquisition stream does).
    if incoming.protocol_state in _FINISHED_PROTOCOL_STATES:
        return _clear_activity_fields(merged)

    return merged


def position_acquisition_active(status: Optional[PositionStatus]) -> bool:
    """Return whether acquisition lifecycle streams indicate an active run."""
    if status is None:
        return False
    state = (status.acquisition_state or "").lower()
    return state in _ACTIVE_ACQUISITION_STATES


def merge_instance_yield(status: PositionStatus, activity: Any) -> PositionStatus:
    """Merge live yield counters from ``stream_instance_activity`` only."""
    yield_summary = getattr(activity, "yield_summary", None)
    if yield_summary is None:
        return status

    updates: dict[str, Any] = {}
    passed = getattr(yield_summary, "basecalled_pass_read_count", None)
    failed = getattr(yield_summary, "basecalled_fail_read_count", None)
    if passed is not None:
        merged_passed = _merge_monotonic_counter(status.passed_reads, passed)
        if merged_passed != status.passed_reads:
            updates["passed_reads"] = merged_passed
    if failed is not None:
        merged_failed = _merge_monotonic_counter(status.failed_reads, failed)
        if merged_failed != status.failed_reads:
            updates["failed_reads"] = merged_failed

    if not updates:
        return status
    return replace(status, **updates)


def merge_instance_activity(
    status: PositionStatus,
    activity: Any,
    *,
    protocol_state_enum: Any = None,
) -> PositionStatus:
    """Backward-compatible merge: yields plus any non-empty run metadata."""
    status = merge_instance_yield(status, activity)
    updates: dict[str, Any] = {}

    flow_cell = getattr(activity, "flow_cell_info", None)
    if flow_cell is not None and getattr(flow_cell, "has_flow_cell", False):
        flow_cell_id = _first_non_empty(
            getattr(flow_cell, "user_specified_flow_cell_id", None),
            getattr(flow_cell, "flow_cell_id", None),
        )
        if flow_cell_id:
            updates["flow_cell_id"] = flow_cell_id
        product_code = _first_non_empty(
            getattr(flow_cell, "user_specified_product_code", None),
            getattr(flow_cell, "product_code", None),
        )
        if product_code:
            updates["flow_cell_product_code"] = product_code

    run = getattr(activity, "protocol_run_info", None)
    if run is not None:
        protocol_run_id = _non_empty_string(getattr(run, "run_id", None))
        if protocol_run_id:
            updates["protocol_run_id"] = protocol_run_id
        protocol_name = _non_empty_string(getattr(run, "protocol_id", None))
        if protocol_name:
            updates["protocol_name"] = protocol_name
        if protocol_state_enum is not None:
            protocol_run_state = _enum_name(
                protocol_state_enum, getattr(run, "state", None)
            )
            if protocol_run_state:
                updates["protocol_run_state"] = protocol_run_state
        run_output = _first_non_empty(
            _non_empty_string(getattr(run, "reported_output_path", None)),
            _non_empty_string(getattr(run, "output_path", None)),
        )
        if run_output:
            updates["output_path"] = run_output

        user_info = getattr(run, "user_info", None)
        if user_info is not None:
            sample_id = _wrapper_value(getattr(user_info, "sample_id", None))
            if sample_id:
                updates["sample_id"] = sample_id
            group_id = _wrapper_value(getattr(user_info, "protocol_group_id", None))
            if group_id:
                updates["protocol_group_id"] = group_id

    if not updates:
        return status
    return replace(status, **updates)


def merge_acquisition_run(
    status: PositionStatus,
    acquisition_run: Any,
    *,
    state_enum: Any,
) -> PositionStatus:
    """Merge fields from ``watch_current_acquisition_run`` into a position snapshot."""
    acquisition_state = _enum_name(state_enum, getattr(acquisition_run, "state", None))
    if not acquisition_state:
        return status

    merged = replace(status, acquisition_state=acquisition_state)
    if acquisition_state in _TERMINAL_ACQUISITION_STATES:
        cleared = _clear_activity_fields(merged)
        return replace(cleared, acquisition_state=acquisition_state)
    return merged


def merge_flow_cell_info(status: PositionStatus, flow_cell: Any) -> PositionStatus:
    """Merge fields from ``get_flow_cell_info`` into a position snapshot."""
    if not getattr(flow_cell, "has_flow_cell", False):
        return status

    updates: dict[str, Any] = {}
    flow_cell_id = _first_non_empty(
        getattr(flow_cell, "user_specified_flow_cell_id", None),
        getattr(flow_cell, "flow_cell_id", None),
    )
    if flow_cell_id:
        updates["flow_cell_id"] = flow_cell_id
    product_code = _first_non_empty(
        getattr(flow_cell, "user_specified_product_code", None),
        getattr(flow_cell, "product_code", None),
    )
    if product_code:
        updates["flow_cell_product_code"] = product_code

    if not updates:
        return status
    return replace(status, **updates)


def merge_protocol_run(
    status: PositionStatus,
    run: Any,
    *,
    protocol_state_enum: Any,
) -> PositionStatus:
    """Merge fields from ``get_current_protocol_run`` into a position snapshot."""
    updates: dict[str, Any] = {}
    protocol_run_id = _non_empty_string(getattr(run, "run_id", None))
    if protocol_run_id:
        updates["protocol_run_id"] = protocol_run_id
    protocol_name = _non_empty_string(getattr(run, "protocol_id", None))
    if protocol_name:
        updates["protocol_name"] = protocol_name
    protocol_run_state = _enum_name(
        protocol_state_enum, getattr(run, "state", None)
    )
    if protocol_run_state:
        updates["protocol_run_state"] = protocol_run_state

    run_output = _first_non_empty(
        _non_empty_string(getattr(run, "reported_output_path", None)),
        _non_empty_string(getattr(run, "output_path", None)),
    )
    if run_output:
        updates["output_path"] = run_output

    user_info = getattr(run, "user_info", None)
    if user_info is not None:
        sample_id = _wrapper_value(getattr(user_info, "sample_id", None))
        if sample_id:
            updates["sample_id"] = sample_id
        group_id = _wrapper_value(getattr(user_info, "protocol_group_id", None))
        if group_id:
            updates["protocol_group_id"] = group_id

    if not updates:
        return status
    return replace(status, **updates)


def merge_output_directories(status: PositionStatus, directories: Any) -> PositionStatus:
    """Attach static output directory paths from ``get_output_directories``."""
    return replace(
        status,
        output_path=_non_empty_string(getattr(directories, "output", None))
        or status.output_path,
        output_reads_path=_non_empty_string(getattr(directories, "reads", None))
        or status.output_reads_path,
        output_logs_path=_non_empty_string(getattr(directories, "log", None))
        or status.output_logs_path,
    )


def simple_protocol_state_name(value: Any) -> str:
    import minknow_api.manager_pb2 as manager_pb2

    return _enum_name(manager_pb2.SimpleProtocolState, value) or "no_protocol_state"


def _wrapper_value(wrapper: Any) -> Optional[str]:
    if wrapper is None:
        return None
    return _non_empty_string(getattr(wrapper, "value", None))


def _non_empty_string(value: Any) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _first_non_empty(*values: Any) -> Optional[str]:
    for value in values:
        text = _non_empty_string(value)
        if text:
            return text
    return None


def _merge_monotonic_counter(
    current: Optional[int],
    new_value: Any,
) -> int:
    """Keep yield counters from decreasing after a yield-stream reconnect."""
    new_int = int(new_value)
    if current is None:
        return new_int
    return max(current, new_int)


def _clear_activity_fields(status: PositionStatus) -> PositionStatus:
    return replace(
        status,
        flow_cell_id=None,
        flow_cell_product_code=None,
        sample_id=None,
        protocol_group_id=None,
        protocol_run_id=None,
        protocol_name=None,
        protocol_run_state=None,
        acquisition_state=None,
        output_path=None,
        output_reads_path=None,
        output_logs_path=None,
        passed_reads=None,
        failed_reads=None,
        connection_error=None,
    )


def _enum_name(enum_type: Any, value: Any) -> Optional[str]:
    if value is None:
        return None
    try:
        name = enum_type.Name(value)
        if name.startswith("STATE_"):
            return name[6:].lower()
        if name.startswith("PROTOCOL_"):
            return name.lower()
        return name.lower()
    except Exception:
        return str(value)
