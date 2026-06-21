"""Shared helpers for mapping MinKNOW protobuf messages to ROBIN models."""

from __future__ import annotations

from dataclasses import replace
from typing import Any, Optional

import minknow_api.manager_pb2 as manager_pb2

from robin.minknow.models import PositionStatus


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


def merge_instance_activity(
    status: PositionStatus,
    activity: Any,
    *,
    protocol_state_enum: Any = None,
) -> PositionStatus:
    """Merge fields from ``stream_instance_activity`` into a position snapshot."""
    updates: dict[str, Any] = {}

    flow_cell = getattr(activity, "flow_cell_info", None)
    if flow_cell is not None and getattr(flow_cell, "has_flow_cell", False):
        updates["flow_cell_id"] = _first_non_empty(
            getattr(flow_cell, "user_specified_flow_cell_id", None),
            getattr(flow_cell, "flow_cell_id", None),
        )
        updates["flow_cell_product_code"] = _first_non_empty(
            getattr(flow_cell, "user_specified_product_code", None),
            getattr(flow_cell, "product_code", None),
        )

    run = getattr(activity, "protocol_run_info", None)
    if run is not None:
        updates["protocol_run_id"] = _non_empty_string(getattr(run, "run_id", None))
        updates["protocol_name"] = _non_empty_string(getattr(run, "protocol_id", None))
        if protocol_state_enum is not None:
            updates["protocol_run_state"] = _enum_name(
                protocol_state_enum, getattr(run, "state", None)
            )
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

    yield_summary = getattr(activity, "yield_summary", None)
    if yield_summary is not None:
        passed = getattr(yield_summary, "basecalled_pass_read_count", None)
        failed = getattr(yield_summary, "basecalled_fail_read_count", None)
        if passed is not None:
            updates["passed_reads"] = int(passed)
        if failed is not None:
            updates["failed_reads"] = int(failed)

    if not updates:
        return status
    return replace(status, **updates)


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
