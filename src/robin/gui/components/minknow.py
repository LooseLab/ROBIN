"""MinKNOW sequencer status panel for the ROBIN GUI."""

from __future__ import annotations

import logging
import queue
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

from nicegui import run, ui

from robin.minknow.config import MinKnowSettings, preset_path_from_environ, workflow_toml_from_environ
from robin.minknow.monitor import (
    MinKnowPollResult,
    fetch_sequencer_status,
    format_poll_summary,
    position_table_rows,
)
from robin.minknow.preset import RobinRunPreset
from robin.minknow.run import (
    MinKnowStartError,
    MinKnowStopError,
    StartRunRequest,
    StopRunRequest,
    start_protocol_run,
    stop_protocol_run,
)
from robin.minknow.sample_id import generate_sample_id_md5
from robin.minknow.stream_monitor import acquire_stream_monitor
from robin.minknow.toml_config import load_minknow_toml
from robin.minknow.workflow_refs import (
    load_workflow_config_for_refs,
    workflow_context_from_runner,
)
from robin.minknow.watch import process_auto_watch, watch_position_run

LOGGER = logging.getLogger(__name__)

_TABLE_COLUMNS = [
    {"name": "position", "label": "Position", "field": "position", "sortable": True, "align": "left"},
    {"name": "state", "label": "State", "field": "state", "sortable": True, "align": "left"},
    {"name": "protocol_state", "label": "Protocol", "field": "protocol_state", "sortable": True, "align": "left"},
    {"name": "sample_id", "label": "Sample ID", "field": "sample_id", "sortable": True, "align": "left"},
    {"name": "protocol_run_id", "label": "Run ID", "field": "protocol_run_id", "sortable": True, "align": "left"},
    {"name": "flow_cell_id", "label": "Flow cell", "field": "flow_cell_id", "sortable": True, "align": "left"},
    {"name": "passed_reads", "label": "Passed reads", "field": "passed_reads", "sortable": True, "align": "left"},
    {"name": "watch_path", "label": "Watch path", "field": "watch_path", "sortable": True, "align": "left"},
    {"name": "actions", "label": "", "field": "actions", "sortable": False, "align": "right"},
]

_ACTIONS_SLOT = """
<q-td key="actions" :props="props">
  <div class="row no-wrap items-center justify-end q-gutter-xs">
    <q-btn
      v-if="props.row.can_watch"
      color="primary"
      size="sm"
      label="Watch"
      icon="folder_open"
      @click="$parent.$emit('watch-run', props.row.position)"
    />
    <q-btn
      v-if="props.row.can_stop"
      color="negative"
      size="sm"
      label="Stop"
      icon="stop"
      @click="$parent.$emit('stop-run', props.row.position)"
    />
    <span v-if="!props.row.can_watch && !props.row.can_stop" class="text-xs text-slate-500">
      No active run
    </span>
  </div>
</q-td>
"""


def add_minknow_sequencer_section(
    *,
    compact: bool = False,
    initial_settings: Optional[MinKnowSettings] = None,
    workflow_runner: Any = None,
    workflow_toml: Optional[Path] = None,
) -> None:
    """Add a live MinKNOW status card driven by manager/activity streams."""
    base_settings = initial_settings or MinKnowSettings.from_environ()
    runner_reference, runner_panel = workflow_context_from_runner(workflow_runner)
    default_workflow_toml = workflow_toml or workflow_toml_from_environ()
    default_preset_path = (
        preset_path_from_environ()
        or default_workflow_toml
        or Path("examples/minknow.example.toml")
    )
    state: dict[str, Any] = {
        "host": base_settings.host,
        "enabled": base_settings.enabled,
        "auto_watch": base_settings.auto_add_paths,
        "release_monitor": None,
        "unsubscribe": None,
        "last_result": MinKnowPollResult(None, None),
        "last_updated": None,
        "preset_path": str(default_preset_path),
        "workflow_toml": str(default_workflow_toml) if default_workflow_toml else "",
        "runner_reference": runner_reference,
        "runner_panel": runner_panel,
        "cached_preset": None,
        "pending_stop": None,
    }
    update_queue: queue.SimpleQueue[MinKnowPollResult] = queue.SimpleQueue()

    def _current_settings() -> MinKnowSettings:
        return MinKnowSettings.from_host(
            host=str(state["host"]).strip() or "localhost",
            enabled=bool(state["enabled"]),
            poll_interval_s=base_settings.poll_interval_s,
            auto_add_paths=bool(state["auto_watch"]),
        )

    def _timestamp() -> str:
        return datetime.now().strftime("%H:%M:%S")

    def _notify(message: str, *, kind: str = "info") -> None:
        try:
            ui.notify(message, type=kind)
        except Exception:
            LOGGER.info("MinKNOW: %s", message)

    def _workflow_config() -> Optional[dict[str, Any]]:
        workflow_path = Path(str(state.get("workflow_toml") or "")).expanduser()
        if workflow_path.is_file():
            return load_workflow_config_for_refs(workflow_path)
        preset_path = Path(str(state.get("preset_path") or "")).expanduser()
        if preset_path.is_file():
            return load_workflow_config_for_refs(preset_path)
        return None

    def _load_preset() -> Optional[RobinRunPreset]:
        path = Path(str(state.get("preset_path") or "")).expanduser()
        if not path.is_file():
            state["cached_preset"] = None
            return None
        try:
            config = load_minknow_toml(
                path,
                workflow_config=_workflow_config(),
                reference=state.get("runner_reference"),
                target_panel=state.get("runner_panel"),
                prefer_workflow=True,
            )
            preset = config.preset
            state["cached_preset"] = preset
            return preset
        except Exception as exc:
            state["cached_preset"] = None
            LOGGER.debug("Failed to load MinKNOW preset", exc_info=True)
            _notify(f"Preset error: {exc}", kind="negative")
            return None

    with ui.element("div").classes("classification-insight-card w-full min-w-0"):
        with ui.column().classes("w-full min-w-0 gap-3 p-2 md:p-3"):
            with ui.row().classes("items-center gap-2 min-w-0 flex-wrap"):
                ui.icon("biotech").classes("classification-insight-icon")
                ui.label("Sequencer (MinKNOW)").classes(
                    "classification-insight-model flex-1 min-w-0"
                )
                stream_indicator = ui.icon("sensors").classes(
                    "classification-insight-icon"
                )
                stream_indicator.set_visibility(False)

            if not compact:
                ui.label(
                    "Live status via MinKNOW streams. Start ROBIN-compliant runs, "
                    "watch BAM output, or stop active protocols."
                ).classes("classification-insight-foot")

            with ui.row().classes("w-full gap-2 flex-wrap items-end"):
                host_input = ui.input(
                    "MinKNOW host",
                    value=state["host"],
                ).props("outlined dense").classes("min-w-[12rem] flex-1")
                enabled_switch = ui.switch(
                    "Monitor",
                    value=state["enabled"],
                ).props("dense")
                auto_watch_switch = ui.switch(
                    "Auto-watch runs",
                    value=state["auto_watch"],
                ).props("dense")
                refresh_button = ui.button(
                    "Refresh now",
                    icon="refresh",
                ).props("flat dense no-caps outline")

            start_controls: dict[str, Any] = {}
            if not compact:
                with ui.expansion("Start ROBIN run", icon="play_arrow").classes("w-full"):
                    with ui.column().classes("w-full gap-2"):
                        preset_input = ui.input(
                            "Preset TOML",
                            value=state["preset_path"],
                        ).props("outlined dense").classes("w-full")
                        preset_input.on(
                            "blur",
                            lambda: state.update(
                                {"preset_path": (preset_input.value or "").strip()}
                            ),
                        )

                        position_input = ui.input(
                            "Position",
                            placeholder="e.g. P2S_000000-A",
                        ).props("outlined dense").classes("w-full")

                        sample_id_input = ui.input(
                            "Sample ID",
                            placeholder="MD5 from Sample ID generator",
                        ).props("outlined dense").classes("w-full")

                        with ui.expansion(
                            "Generate sample ID",
                            icon="fingerprint",
                        ).classes("w-full"):
                            gen_test_id = ui.input("Test ID (required)").props(
                                "outlined dense"
                            ).classes("w-full")
                            gen_first = ui.input("First name").props(
                                "outlined dense"
                            ).classes("w-full")
                            gen_last = ui.input("Last name").props(
                                "outlined dense"
                            ).classes("w-full")
                            gen_dob = ui.input(
                                "Date of birth (YYYY-MM-DD)"
                            ).props("outlined dense").classes("w-full")

                            def _generate_sample_id() -> None:
                                try:
                                    sample_id_input.value = generate_sample_id_md5(
                                        gen_test_id.value or "",
                                        first_name=gen_first.value or "",
                                        last_name=gen_last.value or "",
                                        date_of_birth=gen_dob.value or "",
                                    )
                                    _notify("Sample ID generated.", kind="positive")
                                except ValueError as exc:
                                    _notify(str(exc), kind="warning")

                            ui.button(
                                "Generate",
                                icon="fingerprint",
                                on_click=_generate_sample_id,
                            ).props("flat dense no-caps outline")
                            ui.link(
                                "Open full Sample ID generator",
                                "/sample_id_generator",
                            ).classes("text-xs")

                        experiment_group_input = ui.input(
                            "Experiment group",
                            value="ROBIN_RUN",
                        ).props("outlined dense").classes("w-full")

                        duration_input = ui.number(
                            "Duration (hours)",
                            value=24,
                            min=0.1,
                            step=0.5,
                        ).props("outlined dense").classes("w-full")

                        start_button = ui.button(
                            "Start run…",
                            icon="play_arrow",
                        ).props("color=primary dense no-caps")

                        start_controls.update(
                            {
                                "preset_input": preset_input,
                                "position_input": position_input,
                                "sample_id_input": sample_id_input,
                                "experiment_group_input": experiment_group_input,
                                "duration_input": duration_input,
                                "start_button": start_button,
                            }
                        )

            summary_label = ui.label("Waiting for stream connection…").classes(
                "classification-insight-foot w-full"
            )
            meta_label = ui.label("").classes(
                "text-xs workflow-monitor-meta w-full"
            )
            warning_label = ui.label("").classes(
                "text-xs text-amber-700 dark:text-amber-300 w-full"
            )

            table_container = ui.column().classes("w-full min-w-0")
            with table_container:
                from robin.gui.theme import styled_table

                _, positions_table = styled_table(
                    columns=_TABLE_COLUMNS,
                    rows=[],
                    pagination=5 if compact else 10,
                    class_size="table-xs",
                )
                positions_table.classes("w-full")
                try:
                    positions_table.props(
                        'multi-sort rows-per-page-options="[5,10,20]"'
                    )
                except Exception:
                    pass

                try:
                    positions_table.add_slot("body-cell-actions", _ACTIONS_SLOT)
                except Exception:
                    LOGGER.debug("MinKNOW actions slot failed", exc_info=True)

            error_label = ui.label("").classes(
                "text-sm text-red-600 dark:text-red-400 w-full"
            )

    with ui.dialog() as start_dialog, ui.card().classes("min-w-[20rem]"):
        ui.label("Confirm run start").classes("text-lg font-medium")
        start_confirm_text = ui.label("").classes("text-sm whitespace-pre-wrap")
        with ui.row().classes("w-full justify-end gap-2 mt-2"):
            ui.button("Cancel", on_click=start_dialog.close).props("flat")
            start_confirm_button = ui.button("Start run", icon="play_arrow").props(
                "color=primary"
            )

    with ui.dialog() as stop_dialog, ui.card().classes("min-w-[18rem]"):
        ui.label("Stop sequencing run?").classes("text-lg font-medium")
        stop_confirm_text = ui.label("").classes("text-sm")
        with ui.row().classes("w-full justify-end gap-2 mt-2"):
            ui.button("Cancel", on_click=stop_dialog.close).props("flat")
            stop_confirm_button = ui.button("Stop run", icon="stop").props(
                "color=negative"
            )

    def _position_names() -> list[str]:
        result = state.get("last_result")
        status = result.status if result is not None else None
        if status is None:
            return []
        return [pos.name for pos in status.positions]

    def _find_position_status(name: str):
        result = state.get("last_result")
        status = result.status if result is not None else None
        if status is None:
            return None
        return next((item for item in status.positions if item.name == name), None)

    def _maybe_auto_watch(result: MinKnowPollResult) -> None:
        if not state.get("auto_watch") or result.status is None:
            return
        for position_name, _success, message in process_auto_watch(result.status):
            _notify(f"Auto-watching {position_name}: {message}", kind="positive")

    def _apply_result(result: MinKnowPollResult) -> None:
        state["last_result"] = result
        state["last_updated"] = _timestamp()

        if result.error:
            error_label.set_text(result.error)
            summary_label.set_text(
                format_poll_summary(result, host=state["host"])
            )
            meta_label.set_text("")
            warning_label.set_text("")
            stream_indicator.set_visibility(False)
            positions_table.rows = []
            return

        error_label.set_text("")
        stream_indicator.set_visibility(True)
        status = result.status
        if status is None:
            summary_label.set_text("No MinKNOW data returned.")
            positions_table.rows = []
            return

        summary_label.set_text(
            format_poll_summary(
                result,
                host=f"{status.host}:{status.port}",
                last_updated=state["last_updated"],
            )
        )
        auto_suffix = " · auto-watch on" if state.get("auto_watch") else ""
        meta_label.set_text(
            f"Distribution {status.distribution_version} · "
            f"minknow_api {status.minknow_api_version} · streaming{auto_suffix}"
        )
        warning_label.set_text(status.version_warning or "")
        positions_table.rows = position_table_rows(status)

        if not compact and start_controls:
            pos_input = start_controls.get("position_input")
            if pos_input is not None and not (pos_input.value or "").strip():
                names = _position_names()
                if len(names) == 1:
                    pos_input.value = names[0]

        _maybe_auto_watch(result)

    def _on_stream_update(result: MinKnowPollResult) -> None:
        update_queue.put(result)

    def _drain_stream_updates() -> None:
        latest: Optional[MinKnowPollResult] = None
        while True:
            try:
                latest = update_queue.get_nowait()
            except queue.Empty:
                break
        if latest is not None:
            _apply_result(latest)

    def _detach_stream() -> None:
        unsubscribe = state.get("unsubscribe")
        if unsubscribe is not None:
            try:
                unsubscribe()
            except Exception:
                LOGGER.debug("MinKNOW unsubscribe failed", exc_info=True)
            state["unsubscribe"] = None
        release = state.get("release_monitor")
        if release is not None:
            try:
                release()
            except Exception:
                LOGGER.debug("MinKNOW monitor release failed", exc_info=True)
            state["release_monitor"] = None
        stream_indicator.set_visibility(False)

    def _attach_stream() -> None:
        _detach_stream()
        if not state["enabled"]:
            summary_label.set_text("MinKNOW monitoring is disabled.")
            return
        try:
            monitor, release = acquire_stream_monitor(_current_settings())
        except Exception as exc:
            error_label.set_text(str(exc))
            summary_label.set_text(f"MinKNOW ({state['host']}): not connected")
            return

        state["release_monitor"] = release
        state["unsubscribe"] = monitor.subscribe(_on_stream_update)
        summary_label.set_text(f"Connecting to MinKNOW at {state['host']}…")

    async def _refresh_snapshot() -> None:
        result = await run.io_bound(fetch_sequencer_status, _current_settings())
        _apply_result(result)

    def _on_host_change() -> None:
        state["host"] = (host_input.value or "").strip() or "localhost"

    async def _on_refresh_click() -> None:
        state["host"] = (host_input.value or "").strip() or "localhost"
        state["enabled"] = bool(enabled_switch.value)
        state["auto_watch"] = bool(auto_watch_switch.value)
        await _refresh_snapshot()

    def _on_enabled_change(_event=None) -> None:
        state["enabled"] = bool(enabled_switch.value)
        if state["enabled"]:
            _attach_stream()
        else:
            _detach_stream()
            summary_label.set_text("MinKNOW monitoring is disabled.")

    def _on_auto_watch_change(_event=None) -> None:
        state["auto_watch"] = bool(auto_watch_switch.value)
        if state.get("last_result") and state["last_result"].status is not None:
            meta = meta_label.text or ""
            if state["auto_watch"] and "auto-watch on" not in meta:
                meta_label.set_text(f"{meta} · auto-watch on")
            elif not state["auto_watch"]:
                meta_label.set_text(meta.replace(" · auto-watch on", ""))

    async def _on_host_commit() -> None:
        _on_host_change()
        if state["enabled"]:
            _attach_stream()

    async def _on_watch_run(event) -> None:
        position_name = event.args
        position = _find_position_status(position_name)
        if position is None:
            _notify(f"Position {position_name} not found.", kind="negative")
            return

        success, message = await run.io_bound(watch_position_run, position)
        if success:
            notify_type = "positive"
        elif "already watched" in message.lower():
            notify_type = "warning"
        else:
            notify_type = "negative"
        _notify(message, kind=notify_type)

    def _on_stop_run(event) -> None:
        position_name = event.args
        position = _find_position_status(position_name)
        if position is None:
            _notify(f"Position {position_name} not found.", kind="negative")
            return
        state["pending_stop"] = {
            "position": position_name,
            "run_id": position.protocol_run_id,
            "sample_id": position.sample_id or "—",
        }
        stop_confirm_text.set_text(
            f"Stop the run on {position_name}?\n"
            f"Sample: {position.sample_id or '—'}\n"
            f"Run ID: {position.protocol_run_id or '—'}"
        )
        stop_dialog.open()

    async def _confirm_stop_run() -> None:
        pending = state.get("pending_stop")
        if not pending:
            stop_dialog.close()
            return
        settings = _current_settings()
        try:
            result = await run.io_bound(
                stop_protocol_run,
                settings.auth,
                StopRunRequest(
                    position=pending["position"],
                    protocol_run_id=pending.get("run_id"),
                ),
            )
        except MinKnowStopError as exc:
            _notify(str(exc), kind="negative")
            return
        except Exception as exc:
            _notify(str(exc), kind="negative")
            return
        finally:
            stop_dialog.close()
            state["pending_stop"] = None

        _notify(
            f"Stop requested for {result.position} ({result.protocol_run_id})",
            kind="positive",
        )
        await _refresh_snapshot()

    def _open_start_dialog() -> None:
        if compact or not start_controls:
            return
        state["preset_path"] = (start_controls["preset_input"].value or "").strip()
        preset = _load_preset()
        if preset is None:
            _notify("Preset TOML not found or invalid.", kind="negative")
            return

        position = (start_controls["position_input"].value or "").strip()
        sample_id = (start_controls["sample_id_input"].value or "").strip()
        if not position:
            _notify("Enter a flow cell position.", kind="warning")
            return
        if not sample_id:
            _notify("Enter or generate a sample ID.", kind="warning")
            return

        duration = float(start_controls["duration_input"].value or preset.experiment_duration_hours)
        preset_for_run = preset.with_overrides(experiment_duration_hours=duration)
        errors = preset_for_run.validate()
        if errors:
            _notify("Preset invalid:\n" + "\n".join(errors), kind="negative")
            return

        lines = [
            f"Host: {state['host']}",
            f"Position: {position}",
            f"Sample ID: {sample_id}",
            f"Experiment group: {(start_controls['experiment_group_input'].value or preset.experiment_group).strip()}",
            "",
            *preset_for_run.summary_lines(),
        ]
        state["pending_start"] = {
            "preset": preset_for_run,
            "position": position,
            "sample_id": sample_id,
            "experiment_group": (
                start_controls["experiment_group_input"].value or preset.experiment_group
            ).strip(),
        }
        start_confirm_text.set_text("\n".join(lines))
        start_dialog.open()

    async def _confirm_start_run() -> None:
        pending = state.pop("pending_start", None)
        start_dialog.close()
        if not pending:
            return

        settings = _current_settings()
        request = StartRunRequest(
            preset=pending["preset"],
            position=pending["position"],
            sample_id=pending["sample_id"],
            experiment_group=pending.get("experiment_group"),
        )
        try:
            result = await run.io_bound(
                start_protocol_run,
                settings.auth,
                request,
                validate_models=True,
                check_paths=False,
            )
        except MinKnowStartError as exc:
            _notify(str(exc), kind="negative")
            return
        except Exception as exc:
            _notify(str(exc), kind="negative")
            return

        _notify(
            f"Started run {result.run_id} on {result.position}",
            kind="positive",
        )
        for warning in result.warnings:
            _notify(warning, kind="warning")
        await _refresh_snapshot()

    refresh_button.on_click(_on_refresh_click)
    host_input.on("update:model-value", lambda _e: _on_host_change())
    host_input.on("blur", _on_host_commit)
    enabled_switch.on("update:model-value", _on_enabled_change)
    auto_watch_switch.on("update:model-value", _on_auto_watch_change)
    positions_table.on("watch-run", _on_watch_run)
    positions_table.on("stop-run", _on_stop_run)
    stop_confirm_button.on_click(_confirm_stop_run)
    start_confirm_button.on_click(_confirm_start_run)

    if not compact and start_controls:
        start_controls["start_button"].on_click(_open_start_dialog)

    stream_drain_timer = ui.timer(0.25, _drain_stream_updates, active=True)
    ui.timer(0.2, lambda: _attach_stream() if state["enabled"] else None, once=True)

    def _on_disconnect_cleanup() -> None:
        stream_drain_timer.deactivate()
        _detach_stream()

    try:
        ui.context.client.on_disconnect(_on_disconnect_cleanup)
    except Exception:
        pass
