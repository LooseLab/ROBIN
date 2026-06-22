"""Stream-based MinKNOW monitoring via manager and instance activity APIs."""

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
    merge_instance_activity,
    merge_output_directories,
    merge_position_description,
    position_status_from_description,
)

LOGGER = logging.getLogger(__name__)

Listener = Callable[[MinKnowPollResult], None]
RECONNECT_DELAY_S = 5.0


@dataclass
class _ActivityWorker:
    thread: threading.Thread
    stop_event: threading.Event = field(default_factory=threading.Event)


class MinKnowStreamMonitor:
    """Listen to MinKNOW manager and per-position activity streams."""

    def __init__(self, settings: MinKnowSettings):
        self.settings = settings
        self._listeners: list[Listener] = []
        self._listeners_lock = threading.Lock()
        self._positions: dict[str, PositionStatus] = {}
        self._positions_lock = threading.Lock()
        self._activity_workers: dict[str, _ActivityWorker] = {}
        self._activity_lock = threading.Lock()
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
        self._stop_all_activity_workers()
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
                        self._upsert_position(description, start_activity=True)
                    for description in response.changes:
                        self._upsert_position(description, start_activity=True)
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
                self._stop_all_activity_workers()

            if self._stop_event.is_set():
                break
            time.sleep(RECONNECT_DELAY_S)

    def _upsert_position(self, description: Any, *, start_activity: bool) -> None:
        incoming = position_status_from_description(
            description,
            host=self._host,
            credentials=self._credentials,
        )
        with self._positions_lock:
            current = self._positions.get(incoming.name)
            status = merge_position_description(current, incoming)
            self._positions[incoming.name] = status

        if start_activity and status.running:
            self._ensure_activity_worker(description)
        elif not status.running:
            self._stop_activity_worker(status.name)

    def _remove_position(self, name: str) -> None:
        self._stop_activity_worker(name)
        with self._positions_lock:
            self._positions.pop(name, None)

    def _ensure_activity_worker(self, description: Any) -> None:
        name = description.name
        with self._activity_lock:
            worker = self._activity_workers.get(name)
            if worker is not None and worker.thread.is_alive():
                return
            stop_event = threading.Event()
            thread = threading.Thread(
                target=self._run_activity_stream,
                args=(description, stop_event),
                name=f"minknow-activity-{name}",
                daemon=True,
            )
            self._activity_workers[name] = _ActivityWorker(
                thread=thread, stop_event=stop_event
            )
            thread.start()

    def _stop_activity_worker(self, name: str) -> None:
        with self._activity_lock:
            worker = self._activity_workers.pop(name, None)
        if worker is not None:
            worker.stop_event.set()

    def _stop_all_activity_workers(self) -> None:
        with self._activity_lock:
            workers = list(self._activity_workers.values())
            self._activity_workers.clear()
        for worker in workers:
            worker.stop_event.set()

    def _run_activity_stream(self, description: Any, stop_event: threading.Event) -> None:
        name = description.name
        try:
            from minknow_api.manager import FlowCellPosition

            position = FlowCellPosition(description, self._host, self._credentials)
            connection = position.connect()
        except Exception as exc:
            self._set_connection_error(name, str(exc))
            return

        try:
            import minknow_api.protocol_pb2 as protocol_pb2

            try:
                directories = connection.instance.get_output_directories()
                with self._positions_lock:
                    current = self._positions.get(name)
                    if current is not None:
                        self._positions[name] = merge_output_directories(
                            current, directories
                        )
                self._emit_current()
            except Exception:
                LOGGER.debug(
                    "Could not read output directories for %s", name, exc_info=True
                )

            for activity in connection.instance.stream_instance_activity():
                if stop_event.is_set() or self._stop_event.is_set():
                    break
                with self._positions_lock:
                    current = self._positions.get(name)
                    if current is None:
                        continue
                    self._positions[name] = merge_instance_activity(
                        current,
                        activity,
                        protocol_state_enum=protocol_pb2.ProtocolState,
                    )
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
