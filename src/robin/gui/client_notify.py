"""Send NiceGUI toasts only to clients that still have a Socket.IO connection.

``app.timer`` (used to drain workflow/SNP update queues) has an empty slot stack.
Calling ``ui.notify`` there either raises or can enqueue a message for a client
whose Engine.IO session is already gone after a long SNP run, which surfaces as
an ASGI ``ExceptionGroup`` in Socket.IO ``handle_request``.
"""

from __future__ import annotations

import logging
from typing import Any, Iterable, Optional

logger = logging.getLogger(__name__)

_NOTIFY_TYPES = {"positive", "negative", "warning", "info", "ongoing"}


def iter_connected_clients() -> list[Any]:
    """Return NiceGUI clients that are still connected and not deleted."""
    try:
        from nicegui.client import Client
    except Exception:
        return []
    live: list[Any] = []
    raw_instances = getattr(Client, "instances", None)
    if isinstance(raw_instances, dict):
        values = list(raw_instances.values())
    elif isinstance(raw_instances, (list, tuple, set)):
        values = list(raw_instances)
    else:
        return []
    for client in values:
        try:
            if getattr(client, "is_deleted", False):
                continue
            if not getattr(client, "has_socket_connection", False):
                continue
            live.append(client)
        except Exception:
            continue
    return live


def notify_connected_clients(
    message: str,
    *,
    type: str = "info",
    timeout: Optional[int] = 8000,
    position: str = "top-right",
    clients: Optional[Iterable[Any]] = None,
) -> int:
    """Push ``ui.notify`` to each live client. Returns how many were notified."""
    if not message:
        return 0
    try:
        from nicegui import ui
    except Exception:
        return 0

    notify_type = type if type in _NOTIFY_TYPES else "info"
    sent = 0
    for client in list(clients if clients is not None else iter_connected_clients()):
        try:
            if getattr(client, "is_deleted", False):
                continue
            if not getattr(client, "has_socket_connection", False):
                continue
            with client:
                ui.notify(
                    message,
                    type=notify_type,
                    timeout=timeout,
                    position=position,
                    close_button=True,
                )
            sent += 1
        except Exception:
            logger.debug(
                "Skipping notify for disconnected or invalid NiceGUI client",
                exc_info=True,
            )
    return sent
