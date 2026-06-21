"""Data models for MinKNOW sequencer status."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Optional


@dataclass
class PositionStatus:
    """Status snapshot for a single flow cell position."""

    name: str
    state: str
    protocol_state: str
    running: bool
    device_type: Optional[str] = None
    flow_cell_id: Optional[str] = None
    flow_cell_product_code: Optional[str] = None
    sample_id: Optional[str] = None
    protocol_group_id: Optional[str] = None
    protocol_run_id: Optional[str] = None
    protocol_name: Optional[str] = None
    protocol_run_state: Optional[str] = None
    output_path: Optional[str] = None
    output_reads_path: Optional[str] = None
    output_logs_path: Optional[str] = None
    passed_reads: Optional[int] = None
    failed_reads: Optional[int] = None
    connection_error: Optional[str] = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class SequencerStatus:
    """Status snapshot for a MinKNOW host."""

    host: str
    port: int
    core_version: str
    distribution_version: str
    minknow_api_version: str
    positions: list[PositionStatus] = field(default_factory=list)
    version_warning: Optional[str] = None

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["positions"] = [position.to_dict() for position in self.positions]
        return payload
