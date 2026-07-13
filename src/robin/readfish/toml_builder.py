"""Build readfish experiment TOML from ROBIN presets."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from robin.minknow.preset import RobinRunPreset
from robin.readfish.config import ReadfishConfig, dorado_config_name


def build_readfish_toml_document(
    *,
    preset: RobinRunPreset,
    targets_bed: str,
    minimap2_index: str,
    config: ReadfishConfig,
) -> dict[str, Any]:
    """Return a TOML-serialisable document for ``readfish targets``."""
    if not preset.read_until_filter:
        raise ValueError("read_until_filter is required for readfish adaptive sampling")

    dorado_model = config.dorado_config or dorado_config_name(preset.basecall_simplex_model)
    region = _region_for_filter(
        name="robin_panel",
        filter_mode=preset.read_until_filter,
        targets_bed=str(Path(targets_bed).expanduser()),
        min_chunks=config.min_chunks,
        max_chunks=config.max_chunks,
    )

    mapper_key = "mappy_rs" if config.prom else "mappy"
    mapper_settings: dict[str, Any] = {
        mapper_key: {
            "fn_idx_in": str(Path(minimap2_index).expanduser()),
        }
    }
    if config.prom:
        mapper_settings[mapper_key]["n_threads"] = config.mappy_rs_threads

    return {
        "caller_settings": {
            "dorado": {
                "config": dorado_model,
                "address": config.dorado_address,
            }
        },
        "mapper_settings": mapper_settings,
        "regions": [region],
    }


def _region_for_filter(
    *,
    name: str,
    filter_mode: str,
    targets_bed: str,
    min_chunks: int,
    max_chunks: int,
) -> dict[str, Any]:
    if filter_mode == "enrich":
        on_action = "stop_receiving"
        off_action = "unblock"
    elif filter_mode == "deplete":
        on_action = "unblock"
        off_action = "stop_receiving"
    else:
        raise ValueError(f"Unsupported read_until_filter for readfish: {filter_mode!r}")

    return {
        "name": name,
        "control": False,
        "min_chunks": min_chunks,
        "max_chunks": max_chunks,
        "targets": targets_bed,
        "single_on": on_action,
        "multi_on": on_action,
        "single_off": off_action,
        "multi_off": off_action,
        "no_seq": "proceed",
        "no_map": "proceed",
    }
