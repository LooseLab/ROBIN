"""Ensure readfish TOML artifacts exist when analysing samples with readfish backend."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any, Optional

from robin.readfish.config import DEFAULT_LIVE_TOML_MAX_BIN_WIDTH_BP

LOGGER = logging.getLogger(__name__)

_WORKFLOW_POINTER_NAME = ".robin_workflow_toml"

# Back-compat alias for the default live-update bin-width threshold.
LIVE_TOML_MAX_BIN_WIDTH_BP = DEFAULT_LIVE_TOML_MAX_BIN_WIDTH_BP


def write_workflow_toml_pointer(
    work_dir: str | Path,
    workflow_toml: str | Path,
) -> Path:
    """Persist the active workflow TOML path under ``work_dir`` for Ray workers."""
    work = Path(work_dir).expanduser()
    work.mkdir(parents=True, exist_ok=True)
    pointer = work / _WORKFLOW_POINTER_NAME
    resolved = str(Path(workflow_toml).expanduser().resolve())
    pointer.write_text(resolved + "\n", encoding="utf-8")
    os.environ["ROBIN_WORKFLOW_TOML"] = resolved
    return pointer


def ensure_readfish_toml_for_sample_analysis(
    *,
    sample_id: str,
    work_dir: str | Path,
    master_bed_path: str | Path,
    target_panel: Optional[str] = None,
    reference: Optional[str] = None,
) -> Optional[Path]:
    """Write/update sample ``readfish.toml`` + ``*_live`` when backend is readfish.

    Intended to run after ``master_NNN.bed`` generation during folder analysis so
    adaptive-sampling artifacts appear even without a live MinKNOW start.

    The base ``readfish.toml`` is always written when the backend is readfish. The
    live ``*_live`` targets update is applied only when the sample CNV
    ``bin_width`` is strictly less than ``[readfish].live_toml_max_bin_width_bp``
    (default 1 Mb). Missing ``CNV_dict.npy`` is treated as eligible (no
    coarse-resolution signal yet).
    """
    try:
        context = _resolve_readfish_context(
            work_dir=work_dir,
            target_panel=target_panel,
            reference=reference,
        )
    except Exception:
        LOGGER.debug(
            "Could not resolve readfish workflow context for sample %s",
            sample_id,
            exc_info=True,
        )
        return None

    if context is None:
        return None

    preset, readfish_config = context
    if not preset.readfish_adaptive_sampling_enabled():
        LOGGER.debug(
            "Skipping readfish TOML for sample %s: adaptive_sampling_backend is not readfish",
            sample_id,
        )
        return None

    from robin.readfish.runner import ReadfishStartError, prepare_readfish_toml

    max_bin_width_bp = int(readfish_config.live_toml_max_bin_width_bp)
    allow_live = _bin_width_allows_live_toml_update(
        work_dir=work_dir,
        sample_id=sample_id,
        max_bin_width_bp=max_bin_width_bp,
    )
    if not allow_live:
        LOGGER.info(
            "Skipping readfish live TOML update for sample %s: "
            "CNV bin_width is not strictly below %s bp",
            sample_id,
            max_bin_width_bp,
        )

    try:
        result = prepare_readfish_toml(
            preset=preset,
            config=readfish_config,
            sample_id=sample_id,
            work_directory=Path(work_dir),
            register_live=True,
            master_bed_path=Path(master_bed_path) if allow_live else None,
            work_dir=Path(work_dir) if allow_live else None,
            validate=False,
        )
    except ReadfishStartError as exc:
        LOGGER.warning(
            "Could not write readfish TOML for sample %s during analysis: %s",
            sample_id,
            exc,
        )
        return None
    except Exception:
        LOGGER.warning(
            "Could not write readfish TOML for sample %s during analysis",
            sample_id,
            exc_info=True,
        )
        return None

    live = Path(result.live_toml_path) if result.live_toml_path else None
    LOGGER.info(
        "Analysis wrote readfish TOML for sample %s: %s (live=%s)",
        sample_id,
        result.toml_path,
        live if live and live.is_file() else "pending",
    )
    return live if live and live.is_file() else Path(result.toml_path)


def _sample_cnv_bin_width(
    *,
    work_dir: str | Path,
    sample_id: str,
) -> Optional[int]:
    """Return CNV analysis ``bin_width`` from ``CNV_dict.npy``, or None if unknown."""
    cnv_dict_path = Path(work_dir).expanduser() / sample_id / "CNV_dict.npy"
    if not cnv_dict_path.is_file():
        return None
    try:
        import numpy as np

        cnv_dict = np.load(cnv_dict_path, allow_pickle=True).item()
        bin_width = cnv_dict.get("bin_width")
        if bin_width and int(bin_width) > 0:
            return int(bin_width)
    except Exception:
        LOGGER.debug(
            "Could not load CNV bin_width from %s",
            cnv_dict_path,
            exc_info=True,
        )
    return None


def _bin_width_allows_live_toml_update(
    *,
    work_dir: str | Path,
    sample_id: str,
    max_bin_width_bp: int = LIVE_TOML_MAX_BIN_WIDTH_BP,
) -> bool:
    """True when live TOML targets may be updated for this sample.

    Requires ``bin_width < max_bin_width_bp``. Unknown bin width (no CNV dict yet)
    is allowed so early panel-only master BEDs can still publish live targets.
    """
    bin_width = _sample_cnv_bin_width(work_dir=work_dir, sample_id=sample_id)
    if bin_width is None:
        return True
    return bin_width < int(max_bin_width_bp)


def _resolve_readfish_context(
    *,
    work_dir: str | Path,
    target_panel: Optional[str],
    reference: Optional[str],
) -> Optional[tuple[Any, Any]]:
    from robin.minknow.config import (
        preset_path_from_environ,
        workflow_toml_from_environ,
    )
    from robin.minknow.toml_config import load_minknow_toml
    from robin.readfish.config import ReadfishConfig
    from robin.workflow_config import (
        load_minknow_from_workflow_toml,
        load_workflow_toml,
    )

    for path in _workflow_toml_candidates(work_dir):
        if not path.is_file():
            continue
        try:
            workflow_raw = load_workflow_toml(path)
        except Exception:
            workflow_raw = None

        config = None
        try:
            config = load_minknow_from_workflow_toml(path)
        except Exception:
            config = None
        if config is None:
            try:
                config = load_minknow_toml(
                    path,
                    workflow_config=workflow_raw,
                    prefer_workflow=True,
                )
            except Exception:
                continue

        preset = config.preset
        if preset is None:
            continue

        if target_panel or reference or workflow_raw:
            from robin.minknow.workflow_refs import apply_workflow_refs_to_preset

            preset = apply_workflow_refs_to_preset(
                preset,
                workflow_config=workflow_raw,
                reference=reference,
                target_panel=target_panel,
                prefer_workflow=True,
            )

        readfish = config.readfish if config.readfish is not None else ReadfishConfig()
        return preset, readfish

    preset_path = preset_path_from_environ()
    if preset_path is not None and preset_path.is_file():
        try:
            config = load_minknow_toml(preset_path, prefer_workflow=True)
        except Exception:
            return None
        if config.preset is None:
            return None
        readfish = config.readfish if config.readfish is not None else ReadfishConfig()
        return config.preset, readfish

    _ = workflow_toml_from_environ  # documented discovery path via candidates
    return None


def _workflow_toml_candidates(work_dir: str | Path) -> list[Path]:
    from robin.minknow.config import workflow_toml_from_environ

    candidates: list[Path] = []
    seen: set[str] = set()

    def _add(path: Optional[Path]) -> None:
        if path is None:
            return
        resolved = path.expanduser()
        key = str(resolved)
        if key in seen:
            return
        seen.add(key)
        candidates.append(resolved)

    _add(workflow_toml_from_environ())

    pointer = Path(work_dir).expanduser() / _WORKFLOW_POINTER_NAME
    if pointer.is_file():
        try:
            raw = pointer.read_text(encoding="utf-8").strip()
        except OSError:
            raw = ""
        if raw:
            _add(Path(raw))

    cwd = Path.cwd()
    work = Path(work_dir).expanduser()
    for root in (cwd, work, work.parent):
        for name in (
            "RF.workflow.example.toml",
            "workflow-settings.toml",
            "my_settings.toml",
            "workflow.example.toml",
        ):
            _add(root / name)

    return candidates
