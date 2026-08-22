"""MinKNOW API client for querying sequencer status."""

from __future__ import annotations

import logging
from typing import Any, Optional

import grpc
from packaging.version import Version

from robin.minknow._deps import require_minknow_api
from robin.minknow.auth import MinKnowAuthConfig
from robin.minknow.models import PositionStatus, SequencerStatus
from robin.minknow.parsing import (
    merge_flow_cell_info,
    merge_output_directories,
    merge_protocol_run,
)

LOGGER = logging.getLogger(__name__)


class MinKnowConnectionError(RuntimeError):
    """Raised when MinKNOW cannot be reached or authorisation fails."""


class MinKnowClient:
    """Thin wrapper around ``minknow_api.manager.Manager`` for ROBIN."""

    def __init__(self, auth: MinKnowAuthConfig):
        self.auth = auth
        self._manager = None

    def __enter__(self) -> MinKnowClient:
        self.connect()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def connect(self) -> None:
        require_minknow_api()
        from minknow_api.manager import Manager

        try:
            self._manager = Manager(**self.auth.manager_kwargs())
        except grpc.RpcError as exc:
            raise MinKnowConnectionError(_format_grpc_error(exc)) from exc
        except Exception as exc:
            raise MinKnowConnectionError(str(exc)) from exc

    def close(self) -> None:
        if self._manager is not None:
            try:
                self._manager.close()
            except Exception:
                LOGGER.debug("Error closing MinKNOW manager connection", exc_info=True)
            self._manager = None

    @property
    def manager(self) -> Any:
        if self._manager is None:
            raise RuntimeError("MinKnowClient is not connected; call connect() first.")
        return self._manager

    def get_sequencer_status(self) -> SequencerStatus:
        """Return host version information and per-position status."""
        require_minknow_api()
        import minknow_api

        manager = self.manager
        positions = [
            self._describe_position(position)
            for position in manager.flow_cell_positions()
        ]
        return SequencerStatus(
            host=self.auth.host,
            port=manager.port,
            core_version=manager.core_version,
            distribution_version=manager.version,
            minknow_api_version=minknow_api.__version__,
            positions=positions,
            version_warning=_version_warning(manager, minknow_api.__version__),
        )

    def _describe_position(self, position: Any) -> PositionStatus:
        status = PositionStatus(
            name=position.name,
            state=position.state,
            protocol_state=position.protocol_state,
            running=bool(position.running),
            device_type=getattr(position, "device_type", None),
        )
        if not position.running:
            return status

        try:
            connection = position.connect()
        except Exception as exc:
            status.connection_error = str(exc)
            return status

        self._populate_flow_cell_info(connection, status)
        self._populate_output_directories(connection, status)
        self._populate_current_protocol_run(connection, status)
        return status

    def _populate_flow_cell_info(self, connection: Any, status: PositionStatus) -> None:
        try:
            flow_cell = connection.device.get_flow_cell_info()
        except grpc.RpcError as exc:
            status.connection_error = _format_grpc_error(exc)
            return
        except Exception as exc:
            status.connection_error = str(exc)
            return

        self._apply_status(status, merge_flow_cell_info(status, flow_cell))

    def _populate_output_directories(
        self, connection: Any, status: PositionStatus
    ) -> None:
        try:
            directories = connection.instance.get_output_directories()
        except grpc.RpcError as exc:
            if status.connection_error is None:
                status.connection_error = _format_grpc_error(exc)
            return
        except Exception as exc:
            if status.connection_error is None:
                status.connection_error = str(exc)
            return

        self._apply_status(status, merge_output_directories(status, directories))

    def _populate_current_protocol_run(
        self, connection: Any, status: PositionStatus
    ) -> None:
        try:
            run = connection.protocol.get_current_protocol_run()
        except grpc.RpcError as exc:
            if exc.code() == grpc.StatusCode.FAILED_PRECONDITION:
                return
            if status.connection_error is None:
                status.connection_error = _format_grpc_error(exc)
            return
        except Exception as exc:
            if status.connection_error is None:
                status.connection_error = str(exc)
            return

        self._apply_status(
            status,
            merge_protocol_run(
                status,
                run,
                protocol_state_enum=connection.protocol._pb.ProtocolState,
            ),
        )

    @staticmethod
    def _apply_status(target: PositionStatus, merged: PositionStatus) -> None:
        for field_name in (
            "flow_cell_id",
            "flow_cell_product_code",
            "output_path",
            "output_reads_path",
            "output_logs_path",
            "protocol_run_id",
            "protocol_name",
            "protocol_run_state",
            "sample_id",
            "protocol_group_id",
        ):
            value = getattr(merged, field_name)
            if value is not None:
                setattr(target, field_name, value)


def _version_warning(manager: Any, api_version: str) -> Optional[str]:
    try:
        core_major, core_minor, _core_patch = manager.core_version_components
        pkg = Version(api_version)
    except Exception:
        return None

    if (core_major, core_minor) != (pkg.major, pkg.minor):
        return (
            f"minknow_api {api_version} may be incompatible with "
            f"MinKNOW Core {manager.core_version}; matching minor versions are recommended."
        )
    return None


def _format_grpc_error(exc: grpc.RpcError) -> str:
    code = exc.code().name if exc.code() is not None else "UNKNOWN"
    details = exc.details() or "no details"
    return f"gRPC {code}: {details}"
