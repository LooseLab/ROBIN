"""Stream-based MinKNOW monitoring via manager and per-position gRPC streams."""

from __future__ import annotations

import logging
import threading
import time
from dataclasses import dataclass, field, replace
from typing import Any, Callable, Optional

import grpc

from robin.minknow._deps import require_minknow_api
from robin.minknow.client import _format_grpc_error, _version_warning
from robin.minknow.config import MinKnowSettings
from robin.minknow.models import PositionStatus, SequencerStatus
from robin.minknow.monitor import MinKnowPollResult
from robin.minknow.parsing import (
    merge_acquisition_run,
    merge_flow_cell_info,
    merge_instance_yield,
    merge_position_description,
    merge_protocol_run,
    position_acquisition_active,
    position_status_from_description,
)

LOGGER = logging.getLogger(__name__)

Listener = Callable[[MinKnowPollResult], None]
RECONNECT_DELAY_S = 5.0


@dataclass
class _PositionStreamWorkers:
    threads: list[threading.Thread] = field(default_factory=list)
    stop_event: threading.Event = field(default_factory=threading.Event)


class MinKnowStreamMonitor:
    """Listen to MinKNOW manager and per-position streaming APIs."""

    def __init__(self, settings: MinKnowSettings):
        self.settings = settings
        self._listeners: list[Listener] = []
        self._listeners_lock = threading.Lock()
        self._positions: dict[str, PositionStatus] = {}
        self._positions_lock = threading.Lock()
        self._position_workers: dict[str, _PositionStreamWorkers] = {}
        self._workers_lock = threading.Lock()
        self._host_meta: dict[str, object] = {}
        self._meta_lock = threading.Lock()
        self._stop_event = threading.Event()
        self._thread: Optional[threading.Thread] = None
        self._manager = None
        self._host = settings.auth.host
        self._credentials = None

    def subscribe(self, listener: Listener) -> Callable[[], None]:
        with self._listeners_lock:
            self._listeners.append(listener)
        self._emit_current()

        def unsubscribe() -> None:
            with self._listeners_lock:
                try:
                    self._listeners.remove(listener)
                except ValueError:
                    pass

        return unsubscribe

    def start(self) -> None:
        if self._thread is not None and self._thread.is_alive():
            return
        self._stop_event.clear()
        self._thread = threading.Thread(
            target=self._run_manager_watch,
            name=f"minknow-stream-{self._host}",
            daemon=True,
        )
        self._thread.start()

    def stop(self) -> None:
        self._stop_event.set()
        self._stop_all_position_workers()
        manager = self._manager
        if manager is not None:
            try:
                manager.close()
            except Exception:
                LOGGER.debug("Error closing MinKNOW manager connection", exc_info=True)
        if self._thread is not None:
            self._thread.join(timeout=2.0)
            self._thread = None

    def _run_manager_watch(self) -> None:
        require_minknow_api()
        from minknow_api.manager import Manager

        while not self._stop_event.is_set():
            try:
                manager = Manager(**self.settings.auth.manager_kwargs())
                self._manager = manager
                self._credentials = manager.credentials
                with self._meta_lock:
                    import minknow_api

                    self._host_meta = {
                        "host": self.settings.auth.host,
                        "port": manager.port,
                        "core_version": manager.core_version,
                        "distribution_version": manager.version,
                        "minknow_api_version": minknow_api.__version__,
                        "version_warning": _version_warning(
                            manager, minknow_api.__version__
                        ),
                    }
                self._emit_current()

                for response in manager.rpc.watch_flow_cell_positions():
                    if self._stop_event.is_set():
                        break
                    for description in response.additions:
                        self._upsert_position(description, start_streams=True)
                    for description in response.changes:
                        self._upsert_position(description, start_streams=True)
                    for description in response.removals:
                        self._remove_position(description.name)
                    self._emit_current()
            except grpc.RpcError as exc:
                self._emit_error(_format_grpc_error(exc))
            except Exception as exc:
                LOGGER.debug("MinKNOW stream monitor error", exc_info=True)
                self._emit_error(str(exc))
            finally:
                self._manager = None
                self._stop_all_position_workers()

            if self._stop_event.is_set():
                break
            time.sleep(RECONNECT_DELAY_S)

    def _upsert_position(self, description: Any, *, start_streams: bool) -> None:
        incoming = position_status_from_description(
            description,
            host=self._host,
            credentials=self._credentials,
        )
        with self._positions_lock:
            current = self._positions.get(incoming.name)
            status = merge_position_description(current, incoming)
            self._positions[incoming.name] = status

        if start_streams and status.running:
            self._ensure_position_workers(description)
        elif not status.running and not position_acquisition_active(status):
            self._stop_position_workers(status.name)

    def _remove_position(self, name: str) -> None:
        self._stop_position_workers(name)
        with self._positions_lock:
            self._positions.pop(name, None)

    def _ensure_position_workers(self, description: Any) -> None:
        name = description.name
        with self._workers_lock:
            workers = self._position_workers.get(name)
            if workers is not None and any(
                thread.is_alive() for thread in workers.threads
            ):
                return
            stop_event = threading.Event()
            workers = _PositionStreamWorkers(stop_event=stop_event)
            stream_targets = (
                self._run_protocol_stream,
                self._run_acquisition_stream,
                self._run_yield_stream,
                self._run_flow_cell_stream,
            )
            for target in stream_targets:
                thread = threading.Thread(
                    target=target,
                    args=(description, stop_event),
                    name=f"minknow-{target.__name__}-{name}",
                    daemon=True,
                )
                workers.threads.append(thread)
                thread.start()
            self._position_workers[name] = workers

    def _stop_position_workers(self, name: str) -> None:
        with self._workers_lock:
            workers = self._position_workers.pop(name, None)
        if workers is not None:
            workers.stop_event.set()

    def _stop_all_position_workers(self) -> None:
        with self._workers_lock:
            workers_list = list(self._position_workers.values())
            self._position_workers.clear()
        for workers in workers_list:
            workers.stop_event.set()

    def _connect_position(self, description: Any, name: str) -> Any:
        from minknow_api.manager import FlowCellPosition

        position = FlowCellPosition(description, self._host, self._credentials)
        return position.connect()

    def _update_position(self, name: str, updated: PositionStatus) -> None:
        with self._positions_lock:
            current = self._positions.get(name)
            if current is None:
                return
            self._positions[name] = updated
        self._emit_current()

    def _run_protocol_stream(self, description: Any, stop_event: threading.Event) -> None:
        name = description.name
        try:
            connection = self._connect_position(description, name)
        except Exception as exc:
            self._set_connection_error(name, str(exc))
            return

        try:
            import minknow_api.protocol_pb2 as protocol_pb2

            for run in connection.protocol.watch_current_protocol_run():
                if stop_event.is_set() or self._stop_event.is_set():
                    break
                with self._positions_lock:
                    current = self._positions.get(name)
                    if current is None:
                        continue
                    updated = merge_protocol_run(
                        current,
                        run,
                        protocol_state_enum=protocol_pb2.ProtocolState,
                    )
                    self._positions[name] = updated
                self._emit_current()
        except grpc.RpcError as exc:
            if exc.code() != grpc.StatusCode.CANCELLED:
                self._set_connection_error(name, _format_grpc_error(exc))
        except Exception as exc:
            self._set_connection_error(name, str(exc))

    def _run_acquisition_stream(
        self, description: Any, stop_event: threading.Event
    ) -> None:
        name = description.name
        try:
            connection = self._connect_position(description, name)
        except Exception as exc:
            self._set_connection_error(name, str(exc))
            return

        try:
            import minknow_api.acquisition_pb2 as acquisition_pb2

            for acquisition_run in connection.acquisition.watch_current_acquisition_run():
                if stop_event.is_set() or self._stop_event.is_set():
                    break
                with self._positions_lock:
                    current = self._positions.get(name)
                    if current is None:
                        continue
                    updated = merge_acquisition_run(
                        current,
                        acquisition_run,
                        state_enum=acquisition_pb2.AcquisitionState,
                    )
                    self._positions[name] = updated
                self._emit_current()
        except grpc.RpcError as exc:
            if exc.code() != grpc.StatusCode.CANCELLED:
                self._set_connection_error(name, _format_grpc_error(exc))
        except Exception as exc:
            self._set_connection_error(name, str(exc))

    def _run_yield_stream(self, description: Any, stop_event: threading.Event) -> None:
        name = description.name
        try:
            connection = self._connect_position(description, name)
        except Exception as exc:
            self._set_connection_error(name, str(exc))
            return

        try:
            for activity in connection.instance.stream_instance_activity():
                if stop_event.is_set() or self._stop_event.is_set():
                    break
                with self._positions_lock:
                    current = self._positions.get(name)
                    if current is None:
                        continue
                    updated = merge_instance_yield(current, activity)
                    self._positions[name] = updated
                self._emit_current()
        except grpc.RpcError as exc:
            if exc.code() != grpc.StatusCode.CANCELLED:
                self._set_connection_error(name, _format_grpc_error(exc))
        except Exception as exc:
            self._set_connection_error(name, str(exc))

    def _run_flow_cell_stream(
        self, description: Any, stop_event: threading.Event
    ) -> None:
        name = description.name
        try:
            connection = self._connect_position(description, name)
        except Exception as exc:
            self._set_connection_error(name, str(exc))
            return

        try:
            for flow_cell in connection.device.stream_flow_cell_info():
                if stop_event.is_set() or self._stop_event.is_set():
                    break
                with self._positions_lock:
                    current = self._positions.get(name)
                    if current is None:
                        continue
                    updated = merge_flow_cell_info(current, flow_cell)
                    self._positions[name] = updated
                self._emit_current()
        except grpc.RpcError as exc:
            if exc.code() != grpc.StatusCode.CANCELLED:
                self._set_connection_error(name, _format_grpc_error(exc))
        except Exception as exc:
            self._set_connection_error(name, str(exc))

    def _set_connection_error(self, name: str, message: str) -> None:
        with self._positions_lock:
            current = self._positions.get(name)
            if current is None:
                return
            self._positions[name] = replace(current, connection_error=message)
        self._emit_current()

    def _build_status(self) -> Optional[SequencerStatus]:
        with self._meta_lock:
            if not self._host_meta:
                return None
            with self._positions_lock:
                positions = sorted(self._positions.values(), key=lambda item: item.name)
            return SequencerStatus(
                host=str(self._host_meta["host"]),
                port=int(self._host_meta["port"]),
                core_version=str(self._host_meta["core_version"]),
                distribution_version=str(self._host_meta["distribution_version"]),
                minknow_api_version=str(self._host_meta["minknow_api_version"]),
                positions=positions,
                version_warning=self._host_meta.get("version_warning"),  # type: ignore[arg-type]
            )

    def _emit_current(self) -> None:
        status = self._build_status()
        if status is None:
            return
        self._notify_listeners(MinKnowPollResult(status, None))

    def _emit_error(self, message: str) -> None:
        self._notify_listeners(MinKnowPollResult(None, message))

    def _notify_listeners(self, result: MinKnowPollResult) -> None:
        with self._listeners_lock:
            listeners = list(self._listeners)
        for listener in listeners:
            try:
                listener(result)
            except Exception:
                LOGGER.debug("MinKNOW listener failed", exc_info=True)


# Backward-compatible aliases for tests patching legacy private methods.
MinKnowStreamMonitor._ensure_activity_worker = MinKnowStreamMonitor._ensure_position_workers  # type: ignore[method-assign]
MinKnowStreamMonitor._stop_activity_worker = MinKnowStreamMonitor._stop_position_workers  # type: ignore[method-assign]


@dataclass
class _RegistryEntry:
    monitor: MinKnowStreamMonitor
    refs: int = 0


_registry_lock = threading.Lock()
_registry: dict[str, _RegistryEntry] = {}


def _registry_key(settings: MinKnowSettings) -> str:
    auth = settings.auth
    port = auth.port or 9502
    return f"{auth.host}:{port}"


def acquire_stream_monitor(
    settings: MinKnowSettings,
) -> tuple[MinKnowStreamMonitor, Callable[[], None]]:
    """Return a shared stream monitor for ``settings`` and a release callback."""
    if not settings.enabled:
        raise ValueError("MinKNOW monitoring is disabled.")

    key = _registry_key(settings)
    with _registry_lock:
        entry = _registry.get(key)
        if entry is None:
            monitor = MinKnowStreamMonitor(settings)
            monitor.start()
            entry = _RegistryEntry(monitor=monitor, refs=0)
            _registry[key] = entry
        entry.refs += 1
        monitor = entry.monitor

    def release() -> None:
        with _registry_lock:
            current = _registry.get(key)
            if current is None:
                return
            current.refs -= 1
            if current.refs <= 0:
                current.monitor.stop()
                _registry.pop(key, None)

    return monitor, release
