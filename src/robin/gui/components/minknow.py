"""MinKNOW sequencer status and run control for the ROBIN GUI."""

from __future__ import annotations

import logging
import queue
import sys
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
from robin.minknow.sample_id import (
    build_sample_registration,
    save_sample_registration,
)
from robin.minknow.stream_monitor import acquire_stream_monitor
from robin.minknow.toml_config import MinKnowWorkflowConfig, load_minknow_toml
from robin.minknow.workflow_refs import (
    load_workflow_config_for_refs,
    workflow_context_from_runner,
)
from robin.workflow_config import load_minknow_from_workflow_toml
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


def _workflow_toml_path(
    workflow_toml: Optional[Path],
    *,
    preset_path: Optional[Path] = None,
) -> Optional[Path]:
    if workflow_toml is not None:
        path = workflow_toml.expanduser()
        if path.is_file():
            return path
    if preset_path is not None:
        path = preset_path.expanduser()
        if path.is_file():
            return path
    env_path = workflow_toml_from_environ()
    if env_path is not None and env_path.expanduser().is_file():
        return env_path.expanduser()
    preset_env = preset_path_from_environ()
    if preset_env is not None and preset_env.expanduser().is_file():
        return preset_env.expanduser()
    return None


def _load_workflow_minknow_config(path: Optional[Path]) -> Optional[MinKnowWorkflowConfig]:
    if path is None or not path.is_file():
        return None
    try:
        return load_minknow_from_workflow_toml(path)
    except Exception:
        LOGGER.debug("Failed to load [minknow] from workflow TOML", exc_info=True)
        return None


def add_minknow_sequencer_section(
    *,
    compact: bool = False,
    initial_settings: Optional[MinKnowSettings] = None,
    workflow_runner: Any = None,
    workflow_toml: Optional[Path] = None,
    work_directory: str = "",
) -> None:
    """Add a live MinKNOW status card driven by manager/activity streams."""
    runner_reference, runner_panel = workflow_context_from_runner(workflow_runner)
    resolved_toml = _workflow_toml_path(
        workflow_toml,
        preset_path=preset_path_from_environ(),
    )
    minknow_from_toml = _load_workflow_minknow_config(resolved_toml)
    if minknow_from_toml is not None:
        base_settings = minknow_from_toml.settings
    elif initial_settings is not None:
        base_settings = initial_settings
    else:
        raise RuntimeError(
            "MinKNOW sequencer UI requires workflow TOML or launch-time configuration"
        )

    toml_driven = minknow_from_toml is not None and resolved_toml is not None
    state: dict[str, Any] = {
        "host": base_settings.host,
        "enabled": base_settings.enabled,
        "auto_watch": base_settings.auto_add_paths,
        "release_monitor": None,
        "unsubscribe": None,
        "last_result": MinKnowPollResult(None, None),
        "last_updated": None,
        "workflow_toml": str(resolved_toml) if resolved_toml else "",
        "toml_driven": toml_driven,
        "runner_reference": runner_reference,
        "runner_panel": runner_panel,
        "cached_preset": None,
        "cached_readfish": None,
        "pending_stop": None,
        "selected_position": "",
        "last_position_names": [],
        "position_radio": None,
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
        """Show a GUI toast and mirror the message to the workflow CLI terminal."""
        label = {
            "negative": "ERROR",
            "warning": "WARNING",
            "positive": "OK",
            "info": "INFO",
        }.get(kind, kind.upper())
        text = f"[MinKNOW {label}] {message}"
        stream = sys.stderr if kind in {"negative", "warning"} else sys.stdout
        print(text, file=stream, flush=True)
        if kind == "negative":
            LOGGER.error("%s", message)
        elif kind == "warning":
            LOGGER.warning("%s", message)
        else:
            LOGGER.info("%s", message)
        try:
            ui.notify(message, type=kind)
        except Exception:
            pass

    def _preset_toml_path() -> Optional[Path]:
        workflow_path = Path(str(state.get("workflow_toml") or "")).expanduser()
        if workflow_path.is_file():
            return workflow_path
        return None

    def _workflow_config() -> Optional[dict[str, Any]]:
        workflow_path = _preset_toml_path()
        if workflow_path is not None:
            return load_workflow_config_for_refs(workflow_path)
        return None

    def _load_preset() -> Optional[RobinRunPreset]:
        path = _preset_toml_path()
        if path is None:
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
            state["cached_readfish"] = config.readfish
            if config.settings.host:
                state["host"] = config.settings.host
            state["auto_watch"] = config.settings.auto_add_paths
            return preset
        except Exception as exc:
            state["cached_preset"] = None
            state["cached_readfish"] = None
            LOGGER.debug("Failed to load MinKNOW preset", exc_info=True)
            _notify(f"Preset error: {exc}", kind="negative")
            return None

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

    def _resolve_start_position(preset: RobinRunPreset) -> Optional[str]:
        selected = (state.get("selected_position") or "").strip()
        if selected:
            return selected
        if preset.position:
            return preset.position.strip()
        names = _position_names()
        if len(names) == 1:
            return names[0]
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

            if toml_driven and state.get("workflow_toml"):
                ui.label(
                    f"Defaults from {Path(state['workflow_toml']).name} — "
                    "edit below for this session only (does not change the TOML file)."
                ).classes("text-xs text-slate-500 w-full")

            manual_controls = ui.row().classes("w-full gap-2 flex-wrap items-end")
            with manual_controls:
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

            position_picker_row = None
            position_fallback_input = None
            preset_input = None
            kit_input = None
            simplex_input = None
            modified_input = None
            bam_reads_input = None
            simulation_input = None
            experiment_group_input = None
            duration_input = None
            reference_label = None
            bed_label = None

            toml_preset = _load_preset()
            show_simulation_field = bool(
                toml_preset and toml_preset.simulation_bulk_file
            )

            start_controls: dict[str, Any] = {}

            def _position_has_active_run(name: str) -> bool:
                return False

            def _sync_start_controls_for_position() -> None:
                return None

            def _sync_position_options() -> None:
                return None

            def _apply_selected_position(name: str) -> None:
                name = (name or "").strip()
                if name:
                    state["selected_position"] = name

            if not compact:
                with ui.expansion(
                    "Start ROBIN run",
                    icon="play_arrow",
                    value=True,
                ).classes("w-full"):
                    with ui.column().classes("w-full gap-2"):
                        ui.label(
                            "Choose the flow cell position for this run. Connected "
                            "positions appear as buttons when the monitor is active; "
                            "otherwise type a name or click a row in the table below."
                        ).classes("text-xs text-slate-500 w-full")

                        position_picker_row = ui.row().classes(
                            "w-full gap-2 flex-wrap items-center"
                        )
                        position_fallback_input = ui.input(
                            "Position",
                            placeholder="e.g. P2S_000000-A",
                        ).props("outlined dense").classes("w-full")

                        with ui.expansion(
                            "Run settings",
                            icon="tune",
                            value=True,
                        ).classes("w-full"):
                            with ui.column().classes("w-full gap-2"):
                                if not toml_driven:
                                    ui.label(
                                        "Configure sequencing in workflow TOML "
                                        "([minknow] section) or set host and preset "
                                        "path below."
                                    ).classes("text-xs text-slate-500")
                                    preset_input = ui.input(
                                        "Preset TOML",
                                        value=state.get("workflow_toml") or "",
                                    ).props("outlined dense").classes("w-full")
                                    preset_input.on(
                                        "blur",
                                        lambda: state.update(
                                            {
                                                "workflow_toml": (
                                                    preset_input.value or ""
                                                ).strip()
                                            }
                                        ),
                                    )

                                with ui.row().classes("w-full gap-2 flex-wrap"):
                                    experiment_group_input = ui.input(
                                        "Experiment group",
                                        value="ROBIN_RUN",
                                    ).props("outlined dense readonly").classes(
                                        "flex-1 min-w-[12rem]"
                                    )
                                    duration_input = ui.number(
                                        "Duration (hours)",
                                        value=24,
                                        min=0.1,
                                        step=0.5,
                                    ).props("outlined dense").classes(
                                        "flex-1 min-w-[10rem]"
                                    )

                                kit_input = ui.input(
                                    "Sequencing kit",
                                    value="SQK-LSK114",
                                ).props("outlined dense readonly").classes("w-full")
                                simplex_input = ui.input(
                                    "Basecall simplex model",
                                ).props("outlined dense readonly").classes("w-full")
                                modified_input = ui.input(
                                    "Modified models (comma-separated)",
                                ).props("outlined dense readonly").classes("w-full")
                                with ui.row().classes("w-full gap-2 flex-wrap"):
                                    bam_reads_input = ui.number(
                                        "BAM reads per file",
                                        value=50_000,
                                        min=1,
                                        step=1000,
                                    ).props("outlined dense readonly").classes(
                                        "flex-1 min-w-[12rem]"
                                    )
                                    if show_simulation_field:
                                        simulation_input = ui.input(
                                            "Simulation bulk FAST5",
                                        ).props("outlined dense").classes(
                                            "flex-1 min-w-[12rem]"
                                        )

                                reference_label = ui.label("").classes(
                                    "text-xs text-slate-500 w-full break-all"
                                )
                                bed_label = ui.label("").classes(
                                    "text-xs text-slate-500 w-full break-all"
                                )
                                experiment_group_hint = ui.label("").classes(
                                    "text-xs text-slate-500 w-full"
                                )

                        sample_id_input = ui.input(
                            "Sample ID / MinKNOW RUN ID",
                            placeholder="Registered MD5 or MinKNOW RUN ID",
                        ).props("outlined dense").classes("w-full")

                        with ui.expansion(
                            "Register sample identifiers",
                            icon="fingerprint",
                            value=True,
                        ).classes("w-full"):
                            ui.label(
                                "Enter a MinKNOW RUN ID or generate an MD5 ID, "
                                "optionally with encrypted name/hospital number/notes. "
                                "ROBIN stores the manifest under the work directory and "
                                "links it when this run is detected. Date of birth is "
                                "required when storing encrypted fields."
                            ).classes(
                                "text-xs text-slate-600 dark:text-slate-400 w-full mb-2"
                            )

                            id_mode = ui.toggle(
                                {
                                    "custom": "Use my sample ID",
                                    "md5": "Generate MD5 ID",
                                },
                                value="custom",
                            ).props("no-caps dense").classes("w-full")

                            md5_fields = ui.column().classes("w-full min-w-0 gap-2")
                            with md5_fields:
                                gen_test_id = ui.input(
                                    "Test ID (required for MD5)"
                                ).props("outlined dense").classes("w-full")
                            md5_fields.set_visibility(False)

                            custom_fields = ui.column().classes("w-full min-w-0 gap-2")
                            with custom_fields:
                                custom_run_id = ui.input(
                                    "MinKNOW RUN ID (required)",
                                    placeholder="e.g. HOSP-2024-8841",
                                ).props("outlined dense").classes("w-full font-mono")
                                custom_test_id = ui.input(
                                    "Test ID (optional)"
                                ).props("outlined dense").classes("w-full")

                            gen_first = ui.input("First name (optional)").props(
                                "outlined dense"
                            ).classes("w-full")
                            gen_last = ui.input("Last name (optional)").props(
                                "outlined dense"
                            ).classes("w-full")
                            gen_dob = ui.date_input(
                                "Date of birth (required when encrypting)",
                                value=None,
                            ).classes("w-full")
                            gen_nhs = ui.input(
                                "Hospital number (optional)"
                            ).props("outlined dense").classes("w-full")
                            gen_notes = ui.textarea(
                                "Notes (optional)",
                                placeholder="Free-text notes stored encrypted with identifiers",
                            ).props("outlined dense autogrow").classes("w-full")

                            def _sync_id_mode() -> None:
                                is_md5 = id_mode.value == "md5"
                                md5_fields.set_visibility(is_md5)
                                custom_fields.set_visibility(not is_md5)
                                sample_id_input.value = ""

                            id_mode.on_value_change(lambda _: _sync_id_mode())

                            def _register_identifiers() -> None:
                                try:
                                    registration = build_sample_registration(
                                        mode=id_mode.value or "custom",
                                        custom_sample_id=custom_run_id.value or "",
                                        test_id=(
                                            (gen_test_id.value or "").strip()
                                            if id_mode.value == "md5"
                                            else (custom_test_id.value or "").strip()
                                        ),
                                        first_name=gen_first.value or "",
                                        last_name=gen_last.value or "",
                                        dob=gen_dob.value or "",
                                        nhs_number=gen_nhs.value or "",
                                        notes=gen_notes.value or "",
                                    )
                                except ValueError as exc:
                                    _notify(str(exc), kind="warning")
                                    return

                                sample_id_input.value = registration.sample_id
                                success, msg = save_sample_registration(
                                    work_directory,
                                    registration,
                                )
                                if success:
                                    _notify(
                                        f"Sample ID registered. {msg}",
                                        kind="positive",
                                    )
                                else:
                                    _notify(
                                        f"Sample ID set, but not saved: {msg}",
                                        kind="warning",
                                    )

                            with ui.row().classes("w-full gap-2 flex-wrap"):
                                ui.button(
                                    "Register sample ID",
                                    icon="fingerprint",
                                    on_click=_register_identifiers,
                                ).props("flat dense no-caps outline")
                                ui.link(
                                    "Open Sample ID generator page",
                                    "/sample_id_generator",
                                ).classes("text-xs self-center")

                        start_button = ui.button(
                            "Start run…",
                            icon="play_arrow",
                        ).props("color=primary dense no-caps")
                        start_busy_hint = ui.label("").classes(
                            "text-xs text-amber-700 dark:text-amber-300"
                        )

                        start_controls.update(
                            {
                                "preset_input": preset_input,
                                "sample_id_input": sample_id_input,
                                "id_mode": id_mode,
                                "custom_run_id": custom_run_id,
                                "custom_test_id": custom_test_id,
                                "gen_test_id": gen_test_id,
                                "gen_first": gen_first,
                                "gen_last": gen_last,
                                "gen_dob": gen_dob,
                                "gen_nhs": gen_nhs,
                                "gen_notes": gen_notes,
                                "experiment_group_input": experiment_group_input,
                                "duration_input": duration_input,
                                "kit_input": kit_input,
                                "simplex_input": simplex_input,
                                "modified_input": modified_input,
                                "bam_reads_input": bam_reads_input,
                                "simulation_input": simulation_input,
                                "start_button": start_button,
                                "start_busy_hint": start_busy_hint,
                            }
                        )

                def _apply_selected_position(name: str) -> None:
                    name = (name or "").strip()
                    if not name:
                        return
                    state["selected_position"] = name
                    radio = state.get("position_radio")
                    if radio is not None and radio.value != name:
                        radio.value = name
                    if position_fallback_input is not None:
                        position_fallback_input.value = name
                    _sync_start_controls_for_position()

                def _position_has_active_run(name: str) -> bool:
                    name = (name or "").strip()
                    if not name:
                        return False
                    result = state.get("last_result")
                    status = getattr(result, "status", None) if result is not None else None
                    if status is not None:
                        from robin.minknow.watch import position_has_active_run as _row_active

                        for position in status.positions:
                            if position.name == name:
                                return _row_active(position)
                        return False
                    try:
                        rows = positions_table.rows or []
                    except NameError:
                        return False
                    for row in rows:
                        if row.get("position") == name and row.get("can_stop"):
                            return True
                    return False

                def _sync_start_controls_for_position() -> None:
                    if compact or not start_controls:
                        return
                    button = start_controls.get("start_button")
                    hint = start_controls.get("start_busy_hint")
                    position = (state.get("selected_position") or "").strip()
                    busy = _position_has_active_run(position)
                    if button is not None:
                        if busy:
                            button.disable()
                        else:
                            button.enable()
                    if hint is not None:
                        if busy:
                            hint.set_text(
                                f"{position} already has a run in progress. "
                                "Stop it before starting another."
                            )
                        else:
                            hint.set_text("")

                def _sync_position_options() -> None:
                    if position_picker_row is None:
                        return
                    names = _position_names()
                    preferred = (state.get("selected_position") or "").strip()
                    preset = state.get("cached_preset")
                    if not preferred and preset is not None and preset.position:
                        preferred = preset.position.strip()
                    if not preferred and len(names) == 1:
                        preferred = names[0]

                    if list(names) != state.get("last_position_names"):
                        state["last_position_names"] = list(names)
                        position_picker_row.clear()
                        with position_picker_row:
                            if names:
                                radio_value = (
                                    preferred if preferred in names else names[0]
                                )
                                state["position_radio"] = ui.radio(
                                    names,
                                    value=radio_value,
                                    on_change=lambda e: _apply_selected_position(
                                        e.value
                                    ),
                                ).props("inline").classes("w-full")
                                position_fallback_input.set_visibility(False)
                                _apply_selected_position(radio_value)
                            else:
                                state["position_radio"] = None
                                position_fallback_input.set_visibility(True)
                                if preferred:
                                    position_fallback_input.value = preferred
                                    state["selected_position"] = preferred
                    elif names and state.get("position_radio"):
                        radio = state["position_radio"]
                        if (
                            preferred
                            and preferred in names
                            and radio.value != preferred
                        ):
                            radio.value = preferred

                position_fallback_input.on(
                    "update:model-value",
                    lambda _e: state.update(
                        {
                            "selected_position": (
                                position_fallback_input.value or ""
                            ).strip()
                        }
                    ),
                )

                def _populate_form_from_preset(preset: RobinRunPreset) -> None:
                    experiment_group_input.value = preset.resolve_experiment_group()
                    if preset.append_experiment_group_date_suffix:
                        experiment_group_hint.set_text(
                            "Includes _MON_YY suffix at run start "
                            f"(base name: {preset.experiment_group})"
                        )
                    else:
                        experiment_group_hint.set_text("")
                    duration_input.value = preset.experiment_duration_hours
                    kit_input.value = preset.kit
                    simplex_input.value = preset.basecall_simplex_model
                    modified_input.value = ", ".join(preset.modified_models)
                    bam_reads_input.value = preset.bam_reads_per_file
                    if simulation_input is not None:
                        simulation_input.value = preset.simulation_bulk_file or ""
                    if preset.position:
                        _apply_selected_position(preset.position.strip())
                    if reference_label is not None:
                        reference_label.set_text(
                            "Alignment reference (from workflow): "
                            f"{preset.alignment_reference or '—'}"
                        )
                    if bed_label is not None:
                        bed_label.set_text(
                            f"Stranded panel BED (from workflow): "
                            f"{preset.bed_file or '—'}"
                        )
                    _sync_position_options()

                def _build_preset_from_form(base: RobinRunPreset) -> RobinRunPreset:
                    if simulation_input is not None:
                        simulation_path = (
                            (simulation_input.value or "").strip() or None
                        )
                    else:
                        simulation_path = base.simulation_bulk_file
                    position = (state.get("selected_position") or "").strip()
                    if not position and position_fallback_input is not None:
                        position = (position_fallback_input.value or "").strip()
                    state["selected_position"] = position
                    return base.with_overrides(
                        experiment_duration_hours=float(
                            duration_input.value or base.experiment_duration_hours
                        ),
                        simulation_bulk_file=simulation_path,
                        position=position or None,
                    )

                if toml_preset is not None:
                    _populate_form_from_preset(toml_preset)

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
                    row_key="position",
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

                if not compact:
                    def _on_position_row_click(event) -> None:
                        row = event.args
                        if isinstance(event.args, (list, tuple)) and len(event.args) > 1:
                            row = event.args[1]
                        if not isinstance(row, dict):
                            return
                        name = (row.get("position") or "").strip()
                        if name:
                            _apply_selected_position(name)

                    positions_table.on(
                        "rowClick",
                        _on_position_row_click,
                        [[], ["position"], None],
                    )
                    try:
                        positions_table.props("selection=none")
                    except Exception:
                        pass

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

    with ui.dialog() as stop_dialog, ui.card().classes("min-w-[20rem]"):
        with ui.row().classes("items-center gap-2"):
            ui.icon("warning", color="warning").classes("text-2xl")
            ui.label("Remotely stop sequencing run?").classes("text-lg font-medium")
        stop_confirm_text = ui.label("").classes(
            "text-sm whitespace-pre-wrap text-slate-700 dark:text-slate-300"
        )
        with ui.row().classes("w-full justify-end gap-2 mt-2"):
            ui.button("Cancel", on_click=stop_dialog.close).props("flat")
            stop_confirm_button = ui.button("Stop run", icon="stop").props(
                "color=negative"
            )

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
            if not compact:
                _sync_start_controls_for_position()
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
        if not compact:
            _sync_position_options()
            _sync_start_controls_for_position()

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
            "Are you sure you wish to remotely stop this run? "
            "This will terminate sequencing on the instrument.\n\n"
            f"Position: {position_name}\n"
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

        if start_controls.get("preset_input") is not None:
            state["workflow_toml"] = (
                start_controls["preset_input"].value or ""
            ).strip()

        preset_base = _load_preset()
        if preset_base is None:
            _notify(
                "No [minknow.preset] in workflow TOML. "
                "Add sequencing settings to the file used with robin workflow --toml.",
                kind="negative",
            )
            return

        mode = (start_controls.get("id_mode").value if start_controls.get("id_mode") else None) or "custom"
        sample_id_field = (start_controls["sample_id_input"].value or "").strip()
        custom_run = (
            (start_controls["custom_run_id"].value or "").strip()
            if start_controls.get("custom_run_id") is not None
            else ""
        )
        try:
            if mode == "md5":
                registration = build_sample_registration(
                    mode="md5",
                    test_id=(start_controls["gen_test_id"].value or ""),
                    first_name=(start_controls["gen_first"].value or ""),
                    last_name=(start_controls["gen_last"].value or ""),
                    dob=(start_controls["gen_dob"].value or ""),
                    nhs_number=(start_controls["gen_nhs"].value or ""),
                    notes=(start_controls["gen_notes"].value or ""),
                )
            else:
                registration = build_sample_registration(
                    mode="custom",
                    custom_sample_id=sample_id_field or custom_run,
                    test_id=(start_controls["custom_test_id"].value or ""),
                    first_name=(start_controls["gen_first"].value or ""),
                    last_name=(start_controls["gen_last"].value or ""),
                    dob=(start_controls["gen_dob"].value or ""),
                    nhs_number=(start_controls["gen_nhs"].value or ""),
                    notes=(start_controls["gen_notes"].value or ""),
                )
        except ValueError as exc:
            _notify(str(exc), kind="warning")
            return

        sample_id = registration.sample_id
        start_controls["sample_id_input"].value = sample_id

        ok, msg = save_sample_registration(work_directory, registration)
        if ok:
            _notify(f"Identifiers registered. {msg}", kind="positive")
        else:
            _notify(f"Starting without saved manifest: {msg}", kind="warning")

        preset_for_run = _build_preset_from_form(preset_base)
        position = _resolve_start_position(preset_for_run)
        if not position:
            _notify(
                "Select a flow cell position (or wait for the monitor to list positions).",
                kind="warning",
            )
            return
        if _position_has_active_run(position):
            _notify(
                f"{position} already has a run in progress. "
                "Stop it before starting another.",
                kind="warning",
            )
            return

        experiment_group = preset_for_run.resolve_experiment_group()

        errors = preset_for_run.validate()
        if errors:
            _notify("Preset invalid:\n" + "\n".join(errors), kind="negative")
            return

        lines = [
            f"Host: {state['host']}",
            f"Position: {position}",
            f"Sample ID: {sample_id}",
            f"Experiment group: {experiment_group}",
            "",
            *preset_for_run.summary_lines(),
        ]
        state["pending_start"] = {
            "preset": preset_for_run,
            "position": position,
            "sample_id": sample_id,
            "experiment_group": experiment_group,
            "readfish": state.get("cached_readfish"),
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
            readfish=pending.get("readfish"),
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

        message = f"Started run {result.run_id} on {result.position}"
        if result.readfish_pid is not None:
            message += f" (readfish pid {result.readfish_pid})"
        _notify(
            message,
            kind="positive",
        )
        if result.readfish_pid is not None:
            _notify(
                f"readfish command: {result.readfish_command}",
                kind="info",
            )
            _notify(
                f"readfish dorado: {result.readfish_dorado_config} @ "
                f"{result.readfish_dorado_address}",
                kind="info",
            )
        if result.readfish_log_file:
            _notify(f"readfish log: {result.readfish_log_file}", kind="info")
        if result.readfish_toml_path:
            _notify(f"readfish toml: {result.readfish_toml_path}", kind="info")
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
