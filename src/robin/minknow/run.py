"""Start MinKNOW protocol runs from ROBIN presets."""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import grpc

from robin.minknow._deps import require_minknow_api
from robin.minknow.auth import MinKnowAuthConfig
from robin.minknow.client import MinKnowConnectionError, _format_grpc_error
from robin.minknow.model_resolve import (
    SimplexModelInfo,
    query_basecall_models,
    resolve_preset_simplex_model,
)
from robin.minknow.preset import RobinRunPreset
from robin.minknow.watch import protocol_state_is_inactive
from robin.readfish.config import ReadfishConfig

LOGGER = logging.getLogger(__name__)

# Acquisition must be running before readfish attaches to the device.
_READY_ACQUISITION_STATES = frozenset({"acquisition_running"})
_TERMINAL_ACQUISITION_STATES = frozenset({"acquisition_completed"})


@dataclass(frozen=True)
class BasecallModelsResult:
    """Basecall models available for a flow cell position."""

    product_code: str
    sample_rate: int
    models: list[SimplexModelInfo]


@dataclass(frozen=True)
class StartRunResult:
    """Outcome of starting a MinKNOW protocol run."""

    run_id: str
    position: str
    flow_cell_id: Optional[str]
    protocol_id: str
    sample_id: str
    experiment_group: str
    warnings: tuple[str, ...] = ()
    readfish_pid: Optional[int] = None
    readfish_log_file: Optional[str] = None
    readfish_toml_path: Optional[str] = None
    readfish_command: Optional[str] = None
    readfish_dorado_address: Optional[str] = None
    readfish_dorado_config: Optional[str] = None


@dataclass(frozen=True)
class StartRunRequest:
    """Parameters for starting a single protocol run."""

    preset: RobinRunPreset
    position: str
    sample_id: str
    experiment_group: Optional[str] = None
    readfish: Optional[ReadfishConfig] = None
    work_directory: Optional[str] = None


class MinKnowStartError(RuntimeError):
    """Raised when a protocol run cannot be started."""


class MinKnowStopError(RuntimeError):
    """Raised when a protocol run cannot be stopped."""


def fetch_basecall_models_for_position(
    auth: MinKnowAuthConfig,
    *,
    position: str,
    kit: str,
) -> BasecallModelsResult:
    """List basecall models for the flow cell installed at ``position``."""
    require_minknow_api()
    from minknow_api.manager import Manager
    from minknow_api.tools import protocols

    try:
        manager = Manager(**auth.manager_kwargs())
    except grpc.RpcError as exc:
        raise MinKnowConnectionError(_format_grpc_error(exc)) from exc

    try:
        flow_position = _find_position(manager, position)
        connection = flow_position.connect()
        flow_cell = connection.device.get_flow_cell_info()
        if not getattr(flow_cell, "has_flow_cell", False):
            raise MinKnowStartError(f"No flow cell present in position {position}")

        product_code = getattr(
            flow_cell, "user_specified_product_code", None
        ) or getattr(flow_cell, "product_code", None)
        if not product_code:
            raise MinKnowStartError("Could not determine flow cell product code")

        protocol = protocols.find_protocol(
            connection,
            product_code=product_code,
            kit=kit,
        )
        if protocol is None:
            raise MinKnowStartError(
                f"No matching protocol for kit {kit!r} "
                f"and product code {product_code!r}"
            )

        sample_rate = int(protocol.tags["sample rate"].int_value)
        models, error = query_basecall_models(
            manager,
            product_code=product_code,
            kit=kit,
            sample_rate=sample_rate,
        )
        if error:
            raise MinKnowStartError(error)

        return BasecallModelsResult(
            product_code=product_code,
            sample_rate=sample_rate,
            models=models,
        )
    finally:
        try:
            manager.close()
        except Exception:
            LOGGER.debug("Error closing MinKNOW manager connection", exc_info=True)


def validate_preset_models(
    manager: Any,
    preset: RobinRunPreset,
    *,
    product_code: str,
    sample_rate: int,
) -> list[str]:
    """Validate basecall model names against the connected MinKNOW host."""
    _, _, errors = resolve_preset_simplex_model(
        manager,
        preset,
        product_code=product_code,
        sample_rate=sample_rate,
    )
    return errors


def start_protocol_run(
    auth: MinKnowAuthConfig,
    request: StartRunRequest,
    *,
    validate_models: bool = True,
    check_paths: bool = False,
) -> StartRunResult:
    """Connect to MinKNOW and start a ROBIN-compliant protocol run."""
    preset = request.preset
    experiment_group = preset.resolve_experiment_group(request.experiment_group)
    model_warnings: list[str] = []

    validation_errors = preset.validate(check_paths=check_paths)
    if validation_errors:
        raise MinKnowStartError("; ".join(validation_errors))

    require_minknow_api()
    from minknow_api.manager import Manager
    from minknow_api.tools import protocols

    try:
        manager = Manager(**auth.manager_kwargs())
    except grpc.RpcError as exc:
        raise MinKnowConnectionError(_format_grpc_error(exc)) from exc

    try:
        position = _find_position(manager, request.position)
        connection = position.connect()

        flow_cell = connection.device.get_flow_cell_info()
        if not getattr(flow_cell, "has_flow_cell", False):
            raise MinKnowStartError(
                f"No flow cell present in position {request.position}"
            )

        _ensure_position_idle(connection, request.position)

        product_code = (
            preset.product_code
            or getattr(flow_cell, "user_specified_product_code", None)
            or getattr(flow_cell, "product_code", None)
        )
        if not product_code:
            raise MinKnowStartError("Could not determine flow cell product code")

        protocol = protocols.find_protocol(
            connection,
            product_code=product_code,
            kit=preset.kit,
            config_name=preset.config_name,
        )
        if protocol is None:
            raise MinKnowStartError(
                f"No matching protocol for kit {preset.kit!r} "
                f"and product code {product_code!r}"
            )

        sample_rate = int(protocol.tags["sample rate"].int_value)

        if validate_models:
            preset, model_warnings, model_errors = resolve_preset_simplex_model(
                manager,
                preset,
                product_code=product_code,
                sample_rate=sample_rate,
            )
            for warning in model_warnings:
                LOGGER.warning("MinKNOW model resolution: %s", warning)
            if model_errors:
                raise MinKnowStartError("; ".join(model_errors))

        if preset.product_code:
            connection.device.set_user_specified_product_code(code=preset.product_code)

        basecalling_args = None
        if preset.enable_basecalling:
            alignment_args = protocols.AlignmentArgs(
                reference_files=[preset.alignment_reference],
                bed_file=preset.bed_file,
            )
            basecalling_args = protocols.BasecallingArgs(
                simplex_model=preset.basecall_simplex_model,
                modified_models=list(preset.modified_models) or None,
                stereo_model=None,
                barcoding=None,
                alignment=alignment_args,
                min_qscore=_default_min_qscore(
                    manager, preset, sample_rate, product_code
                ),
            )

        read_until_args = None
        if preset.minknow_adaptive_sampling_enabled():
            read_until_args = protocols.ReadUntilArgs(
                filter_type=preset.read_until_filter,
                reference_files=[preset.effective_read_until_reference()],
                bed_file=preset.effective_read_until_bed_file(),
                first_channel=None,
                last_channel=None,
            )

        bam_arguments = None
        if preset.enable_bam:
            bam_arguments = protocols.OutputArgs(
                reads_per_file=preset.bam_reads_per_file,
                batch_duration=preset.bam_batch_duration,
            )

        stop_criteria = protocols.CriteriaValues(
            runtime=int(preset.experiment_duration_hours * 60 * 60)
        )

        start_kwargs: dict[str, Any] = {
            "disable_active_channel_selection": False,
            "mux_scan_period": preset.mux_scan_period,
            "stop_criteria": stop_criteria,
        }
        start_kwargs.update(_simulation_start_kwargs(preset))

        run_id = protocols.start_protocol(
            connection,
            identifier=protocol.identifier,
            sample_id=request.sample_id,
            experiment_group=experiment_group,
            barcode_info=None,
            basecalling=basecalling_args,
            read_until=read_until_args,
            fastq_arguments=None,
            fast5_arguments=None,
            pod5_arguments=None,
            bam_arguments=bam_arguments,
            **start_kwargs,
        )

        readfish_pid: Optional[int] = None
        readfish_log_file: Optional[str] = None
        readfish_toml_path: Optional[str] = None
        readfish_command: Optional[str] = None
        readfish_dorado_address: Optional[str] = None
        readfish_dorado_config: Optional[str] = None
        if preset.readfish_adaptive_sampling_enabled():
            from robin.readfish.config import ReadfishConfig
            from robin.readfish.runner import ReadfishStartError, start_readfish_targets

            readfish_config = (
                request.readfish if request.readfish is not None else ReadfishConfig()
            )
            print(
                "[readfish] MinKNOW protocol started "
                f"(run_id={run_id}); waiting for acquisition to run…",
                flush=True,
            )
            wait_for_acquisition_running(
                connection,
                protocol_run_id=run_id,
                position_name=position.name,
                timeout_seconds=readfish_config.start_wait_timeout_seconds,
                poll_seconds=readfish_config.start_wait_poll_seconds,
            )
            print(
                "[readfish] Acquisition running; launching readfish adaptive sampling…",
                flush=True,
            )
            try:
                readfish_result = start_readfish_targets(
                    preset=preset,
                    config=readfish_config,
                    auth=auth,
                    position=position.name,
                    sample_id=request.sample_id,
                    experiment_group=experiment_group,
                    work_directory=(
                        Path(request.work_directory) if request.work_directory else None
                    ),
                )
            except ReadfishStartError as exc:
                raise MinKnowStartError(
                    f"MinKNOW protocol started (run_id={run_id}) but readfish failed: {exc}"
                ) from exc
            readfish_pid = readfish_result.pid
            readfish_log_file = readfish_result.log_file
            readfish_toml_path = readfish_result.toml_path
            readfish_command = " ".join(readfish_result.command)
            readfish_dorado_address = readfish_result.dorado_address
            readfish_dorado_config = readfish_result.dorado_config
        elif preset.minknow_adaptive_sampling_enabled():
            print(
                "[readfish] Adaptive sampling backend is minknow (native Read Until); "
                "readfish will not be launched.",
                flush=True,
            )
        else:
            print(
                "[readfish] Adaptive sampling is off; readfish will not be launched.",
                flush=True,
            )

        return StartRunResult(
            run_id=run_id,
            position=position.name,
            flow_cell_id=getattr(flow_cell, "flow_cell_id", None),
            protocol_id=protocol.identifier,
            sample_id=request.sample_id,
            experiment_group=experiment_group,
            warnings=tuple(model_warnings),
            readfish_pid=readfish_pid,
            readfish_log_file=readfish_log_file,
            readfish_toml_path=readfish_toml_path,
            readfish_command=readfish_command,
            readfish_dorado_address=readfish_dorado_address,
            readfish_dorado_config=readfish_dorado_config,
        )
    finally:
        try:
            manager.close()
        except Exception:
            LOGGER.debug("Error closing MinKNOW manager", exc_info=True)


@dataclass(frozen=True)
class StopRunRequest:
    """Parameters for stopping a protocol run."""

    position: str
    protocol_run_id: Optional[str] = None


@dataclass(frozen=True)
class StopRunResult:
    """Outcome of stopping a protocol run."""

    position: str
    protocol_run_id: str
    protocol_state: Optional[str] = None
    waited: bool = False


def stop_protocol_run(
    auth: MinKnowAuthConfig,
    request: StopRunRequest,
    *,
    wait: bool = False,
) -> StopRunResult:
    """Stop a running protocol on a MinKNOW position."""
    require_minknow_api()
    from minknow_api.manager import Manager

    try:
        manager = Manager(**auth.manager_kwargs())
    except grpc.RpcError as exc:
        raise MinKnowConnectionError(_format_grpc_error(exc)) from exc

    protocol_state: Optional[str] = None
    try:
        position = _find_position(manager, request.position, error_cls=MinKnowStopError)
        connection = position.connect()
        protocol_run_id = (request.protocol_run_id or "").strip()

        if not protocol_run_id:
            try:
                run = connection.protocol.get_current_protocol_run()
            except grpc.RpcError as exc:
                if exc.code() == grpc.StatusCode.FAILED_PRECONDITION:
                    raise MinKnowStopError(
                        f"No active protocol run on position {request.position}"
                    ) from exc
                raise MinKnowStopError(_format_grpc_error(exc)) from exc
            protocol_run_id = getattr(run, "run_id", None) or ""
            if not protocol_run_id:
                raise MinKnowStopError(
                    f"No protocol run ID reported for position {request.position}"
                )

        connection.protocol.stop_protocol(protocol_run_id=protocol_run_id)

        if wait:
            run_info = connection.protocol.wait_for_finished(run_id=protocol_run_id)
            state_enum = connection.protocol._pb.ProtocolState
            protocol_state = _enum_name(state_enum, getattr(run_info, "state", None))
    finally:
        try:
            manager.close()
        except Exception:
            LOGGER.debug("Error closing MinKNOW manager", exc_info=True)

    return StopRunResult(
        position=request.position,
        protocol_run_id=protocol_run_id,
        protocol_state=protocol_state,
        waited=wait,
    )


def wait_for_acquisition_running(
    connection: Any,
    *,
    protocol_run_id: str,
    position_name: str,
    timeout_seconds: float = 600.0,
    poll_seconds: float = 10.0,
) -> str:
    """Block until MinKNOW reports ``ACQUISITION_RUNNING`` for the position.

    readfish must attach after acquisition has actually started; launching during
    temperature wait / mux scan / protocol setup causes connection failures.
    """
    timeout = max(1.0, float(timeout_seconds))
    poll = max(0.2, float(poll_seconds))
    deadline = time.monotonic() + timeout
    last_acq_state: Optional[str] = None
    last_protocol_state: Optional[str] = None

    while True:
        protocol_state = _current_protocol_state(connection)
        last_protocol_state = protocol_state or last_protocol_state
        if protocol_state_is_inactive(protocol_state):
            raise MinKnowStartError(
                f"Protocol on {position_name} ended before acquisition started "
                f"(run_id={protocol_run_id}, protocol_state={protocol_state})."
            )

        acq_state = _current_acquisition_state(connection)
        if acq_state:
            last_acq_state = acq_state
        if acq_state in _READY_ACQUISITION_STATES:
            print(
                f"[readfish] Acquisition ready on {position_name} "
                f"(state={acq_state}).",
                flush=True,
            )
            return acq_state
        if acq_state in _TERMINAL_ACQUISITION_STATES:
            raise MinKnowStartError(
                f"Acquisition on {position_name} finished before readfish could start "
                f"(run_id={protocol_run_id}, acquisition_state={acq_state})."
            )

        if time.monotonic() >= deadline:
            detail = f"run_id={protocol_run_id}"
            if last_protocol_state:
                detail += f", protocol_state={last_protocol_state}"
            if last_acq_state:
                detail += f", acquisition_state={last_acq_state}"
            else:
                detail += ", acquisition_state=none"
            raise MinKnowStartError(
                f"Timed out after {timeout:.0f}s waiting for acquisition on "
                f"{position_name} ({detail})."
            )

        print(
            "[readfish] Waiting for acquisition…"
            f" protocol={last_protocol_state or 'unknown'}"
            f" acquisition={last_acq_state or 'none'}",
            flush=True,
        )
        time.sleep(poll)


def _current_protocol_state(connection: Any) -> Optional[str]:
    try:
        run = connection.protocol.get_current_protocol_run()
    except grpc.RpcError as exc:
        if exc.code() == grpc.StatusCode.FAILED_PRECONDITION:
            return None
        raise MinKnowStartError(_format_grpc_error(exc)) from exc
    except Exception as exc:
        raise MinKnowStartError(
            f"Could not read current protocol state: {exc}"
        ) from exc

    return _enum_name(
        getattr(getattr(connection.protocol, "_pb", None), "ProtocolState", None),
        getattr(run, "state", None),
    )


def _current_acquisition_state(connection: Any) -> Optional[str]:
    try:
        run = connection.acquisition.get_current_acquisition_run()
    except grpc.RpcError as exc:
        if exc.code() == grpc.StatusCode.FAILED_PRECONDITION:
            return None
        raise MinKnowStartError(_format_grpc_error(exc)) from exc
    except Exception as exc:
        raise MinKnowStartError(
            f"Could not read current acquisition state: {exc}"
        ) from exc

    state_enum = getattr(
        getattr(connection.acquisition, "_pb", None), "AcquisitionState", None
    )
    if state_enum is None:
        try:
            import minknow_api.acquisition_pb2 as acquisition_pb2

            state_enum = acquisition_pb2.AcquisitionState
        except Exception:
            state_enum = None
    return _enum_name(state_enum, getattr(run, "state", None))


def _enum_name(enum_type: Any, value: Any) -> Optional[str]:
    if value is None:
        return None
    try:
        return enum_type.Name(value).lower()
    except Exception:
        return str(value)


def _ensure_position_idle(connection: Any, position_name: str) -> None:
    """Refuse to start when the position already has a protocol run in progress."""
    try:
        run = connection.protocol.get_current_protocol_run()
    except grpc.RpcError as exc:
        if exc.code() == grpc.StatusCode.FAILED_PRECONDITION:
            return
        raise MinKnowStartError(_format_grpc_error(exc)) from exc
    except Exception as exc:
        raise MinKnowStartError(
            f"Could not check for an active run on {position_name}: {exc}"
        ) from exc

    run_id = getattr(run, "run_id", None)
    run_id_text = str(run_id).strip() if run_id is not None else ""
    if not run_id_text:
        return

    state = _enum_name(
        getattr(getattr(connection.protocol, "_pb", None), "ProtocolState", None),
        getattr(run, "state", None),
    )
    if protocol_state_is_inactive(state):
        return

    sample_id = None
    user_info = getattr(run, "user_info", None)
    if user_info is not None:
        sample_wrapper = getattr(user_info, "sample_id", None)
        sample_id = getattr(sample_wrapper, "value", None) or sample_wrapper
        if sample_id is not None:
            sample_id = str(sample_id).strip() or None

    detail = f"run_id={run_id_text}"
    if sample_id:
        detail += f", sample_id={sample_id}"
    if state:
        detail += f", protocol_state={state}"
    raise MinKnowStartError(
        f"Position {position_name} already has a run in progress ({detail}). "
        "Stop the current run before starting another."
    )


def _simulation_start_kwargs(preset: RobinRunPreset) -> dict[str, Any]:
    """Build kwargs for MinKNOW simulated bulk FAST5 playback."""
    if not preset.simulation_bulk_file:
        return {}

    bulk = Path(preset.simulation_bulk_file).expanduser()
    # minknow_api checks simulation_path.exists() on the ROBIN host; when MinKNOW
    # is remote, pass --simulation via args so only the sequencer needs the file.
    if bulk.is_file():
        return {"simulation_path": bulk}
    return {"args": ["--simulation", str(bulk)]}


def _find_position(
    manager: Any,
    name: str,
    *,
    error_cls: type[RuntimeError] = MinKnowStartError,
) -> Any:
    for position in manager.flow_cell_positions():
        if position.name == name:
            return position
    available = ", ".join(pos.name for pos in manager.flow_cell_positions())
    raise error_cls(f"Position {name!r} not found. Available: {available or 'none'}")


def _default_min_qscore(
    manager: Any,
    preset: RobinRunPreset,
    sample_rate: int,
    product_code: str,
) -> Optional[int]:
    from minknow_api.tools import protocols

    try:
        configs = manager.find_basecall_configurations(
            product_code, preset.kit, sample_rate
        )
        simplex = protocols.find_simplex_model(configs, preset.basecall_simplex_model)
        return int(simplex.default_q_score_cutoff)
    except Exception:
        LOGGER.debug("Could not resolve default q-score cutoff", exc_info=True)
        return None


def _available_simplex_models(configs: Any) -> list[str]:
    from robin.minknow.model_resolve import available_simplex_models

    return available_simplex_models(configs)
