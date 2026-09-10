from __future__ import annotations

from robin.gui.client_notify import (
    iter_connected_clients,
    notify_connected_clients,
    run_javascript_when_connected,
    schedule_after_page_sent,
)


class _FakeClient:
    def __init__(self, *, deleted: bool = False, connected: bool = True) -> None:
        self.is_deleted = deleted
        self.has_socket_connection = connected
        self.entered = False
        self.connect_handlers = []

    def __enter__(self):
        self.entered = True
        return self

    def __exit__(self, *args):
        return False

    def on_connect(self, handler):
        self.connect_handlers.append(handler)


def test_iter_connected_clients_filters_dead(monkeypatch) -> None:
    live = _FakeClient(deleted=False, connected=True)
    dead = _FakeClient(deleted=True, connected=True)
    stale = _FakeClient(deleted=False, connected=False)

    from nicegui.client import Client

    monkeypatch.setattr(
        Client, "instances", {"a": live, "b": dead, "c": stale}, raising=True
    )
    assert iter_connected_clients() == [live]


def test_notify_connected_clients_skips_disconnected(monkeypatch) -> None:
    live = _FakeClient(connected=True)
    stale = _FakeClient(connected=False)
    notified: list[str] = []

    class _Ui:
        @staticmethod
        def notify(message, **kwargs):
            notified.append(message)

    monkeypatch.setattr("nicegui.ui", _Ui)

    sent = notify_connected_clients("SNP batch finished", clients=[live, stale])
    assert sent == 1
    assert notified == ["SNP batch finished"]
    assert live.entered is True
    assert stale.entered is False


def test_notify_connected_clients_skips_empty_message() -> None:
    assert notify_connected_clients("") == 0


def test_run_javascript_when_connected_waits_for_socket(monkeypatch) -> None:
    client = _FakeClient(connected=False)
    ran = []

    class _Context:
        client = None

    class _Ui:
        context = _Context()

        @staticmethod
        def run_javascript(code, **kwargs):
            ran.append(code)

    _Ui.context.client = client
    monkeypatch.setattr("nicegui.ui", _Ui)

    run_javascript_when_connected("window.ready = true")
    assert ran == []
    assert client.connect_handlers
    client.connect_handlers[0]()
    assert ran == ["window.ready = true"]


def test_run_javascript_when_connected_runs_immediately_if_connected(monkeypatch) -> None:
    client = _FakeClient(connected=True)
    ran = []

    class _Context:
        client = None

    class _Ui:
        context = _Context()

        @staticmethod
        def run_javascript(code, **kwargs):
            ran.append(code)

    _Ui.context.client = client
    monkeypatch.setattr("nicegui.ui", _Ui)

    run_javascript_when_connected("window.ready = true")
    assert ran == ["window.ready = true"]
    assert client.connect_handlers == []


def test_schedule_after_page_sent_uses_timer(monkeypatch) -> None:
    scheduled = []

    class _Ui:
        @staticmethod
        def timer(delay, callback, once=False):
            scheduled.append((delay, once))
            callback()

    monkeypatch.setattr("nicegui.ui", _Ui)
    hits = []
    schedule_after_page_sent(lambda: hits.append("ok"), delay_s=0.05)
    assert scheduled == [(0.05, True)]
    assert hits == ["ok"]
