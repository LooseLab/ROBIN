"""Workflow integration for ITD / insertion hotspot calling."""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

from robin.analysis.itd_work import (
    accumulate_itd_candidates,
    process_bam_itd_with_staging,
    resolve_itd_hotspots,
)
from robin.logging_config import get_job_logger

logger = logging.getLogger(__name__)

ITD_CONFIG_FILENAME = "_robin_itd_config.json"


def _get_itd_batch_size(default: int = 20) -> int:
    try:
        from robin.workflow_ray import BATCH_CONFIG

        return int(BATCH_CONFIG.get("itd", {}).get("max_batch_size", default))
    except Exception:
        return default


# Optional defaults from workflow TOML ``[itd]`` section (same-process only).
_ITD_TOML_DEFAULTS: Dict[str, Any] = {}


def configure_itd_defaults(
    config: Optional[Dict[str, Any]],
    *,
    work_dir: Optional[str | Path] = None,
) -> None:
    """
    Apply ``[itd]`` settings from workflow TOML.

    Also writes ``{work_dir}/_robin_itd_config.json`` when ``work_dir`` is set so
    Ray workers (separate processes) can load the same settings.
    """
    global _ITD_TOML_DEFAULTS
    _ITD_TOML_DEFAULTS = dict(config or {})
    if work_dir is None or not _ITD_TOML_DEFAULTS:
        return
    try:
        out_dir = Path(work_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / ITD_CONFIG_FILENAME
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(_ITD_TOML_DEFAULTS, handle, indent=2)
        logger.info("Wrote ITD config for workers: %s (%s)", path, _ITD_TOML_DEFAULTS)
    except Exception as exc:
        logger.warning("Could not write ITD config to work_dir: %s", exc)


def _load_itd_config_file(work_dir: Optional[str | Path]) -> Dict[str, Any]:
    if not work_dir:
        return {}
    path = Path(work_dir) / ITD_CONFIG_FILENAME
    if not path.is_file():
        return {}
    try:
        with open(path, encoding="utf-8") as handle:
            raw = json.load(handle)
        return raw if isinstance(raw, dict) else {}
    except Exception as exc:
        logger.debug("Could not read ITD config %s: %s", path, exc)
        return {}


def _load_active_hotspots(
    target_panel: str,
    metadata: Dict[str, Any],
    *,
    work_dir: Optional[str | Path] = None,
):
    nested = metadata.get("itd") if isinstance(metadata.get("itd"), dict) else {}
    # Priority: job metadata > in-process TOML defaults > work_dir config file.
    cfg = {
        **_load_itd_config_file(work_dir),
        **_ITD_TOML_DEFAULTS,
        **nested,
    }

    region_mode = (
        metadata.get("itd_region_mode")
        or cfg.get("region_mode")
        or cfg.get("itd_region_mode")
        or "hotspots"
    )
    hotspots_path = (
        metadata.get("itd_hotspots_path")
        or cfg.get("hotspots_path")
        or cfg.get("itd_hotspots_path")
    )
    hotspots = resolve_itd_hotspots(
        target_panel,
        region_mode=str(region_mode).strip().lower(),
        hotspots_path=hotspots_path,
        panel_min_length=int(
            metadata.get("itd_panel_min_length", cfg.get("panel_min_length", 4))
        ),
        panel_min_frequency=float(
            metadata.get(
                "itd_panel_min_frequency", cfg.get("panel_min_frequency", 0.05)
            )
        ),
        panel_min_supporting_reads=int(
            metadata.get(
                "itd_panel_min_supporting_reads",
                cfg.get("panel_min_supporting_reads", 3),
            )
        ),
    )
    logger.info(
        "ITD region_mode=%s → %d scan windows for panel %s",
        str(region_mode).strip().lower(),
        len(hotspots),
        target_panel,
    )
    return hotspots


def _resolve_sample_id(job) -> str:
    try:
        sample_id = job.context.get_sample_id()
        if sample_id and sample_id != "unknown":
            return sample_id
    except Exception:
        pass
    bam_md = job.context.metadata.get("bam_metadata", {}) or {}
    return str(bam_md.get("sample_id") or "unknown")


def process_multiple_bams_itd(
    bam_paths: List[str],
    *,
    work_dir: str,
    sample_id: str,
    target_panel: str,
    hotspots,
    batch_size: int,
    log: logging.Logger,
) -> Dict[str, Any]:
    """Stage each BAM then force-accumulate ITD candidates."""
    staged = 0
    for bam_path in bam_paths:
        try:
            stats, _ = process_bam_itd_with_staging(
                bam_path,
                work_dir,
                sample_id,
                hotspots,
                batch_size=batch_size,
            )
            staged += 1
            log.info(
                "ITD staged %s (%d count rows, pending=%s)",
                os.path.basename(bam_path),
                stats.get("rows", 0),
                stats.get("pending"),
            )
        except Exception as exc:
            log.exception("ITD staging failed for %s: %s", bam_path, exc)

    accumulation = accumulate_itd_candidates(
        work_dir,
        sample_id,
        hotspots,
        force=True,
        batch_size=batch_size,
    )
    return {
        "files_processed": staged,
        "total_files": len(bam_paths),
        "accumulation": accumulation,
        "events": accumulation.get("events", 0),
        "events_path": accumulation.get("events_path"),
        "summary_path": accumulation.get("summary_path"),
        "genes": sorted(hotspots),
    }


def itd_handler(job, work_dir=None, target_panel=None):
    """
    Handler for ITD / insertion analysis jobs.

    Args:
        job: Workflow Job
        work_dir: Output root (sample dirs created underneath)
        target_panel: Active gene panel name (required)
    """
    if not target_panel:
        raise ValueError("target_panel is required for ITD analysis")

    log = get_job_logger(str(job.job_id), "itd", job.context.filepath)
    batch_size = _get_itd_batch_size()

    if work_dir is None:
        work_dir = os.path.dirname(job.context.filepath) or "."
    os.makedirs(work_dir, exist_ok=True)

    try:
        hotspots = _load_active_hotspots(
            target_panel, job.context.metadata, work_dir=work_dir
        )
    except Exception as exc:
        log.error("Failed to load ITD hotspots for panel %s: %s", target_panel, exc)
        job.context.add_error("itd_analysis", str(exc))
        return

    if not hotspots:
        log.warning(
            "No ITD hotspots overlap panel %s; writing empty results", target_panel
        )

    batched_job = job.context.metadata.get("_batched_job")
    if batched_job:
        sample_id = batched_job.get_sample_id()
        filepaths = batched_job.get_filepaths()
        log.info(
            "ITD batch: %d files for sample %s (panel=%s, %d genes)",
            len(filepaths),
            sample_id,
            target_panel,
            len(hotspots),
        )
        result = process_multiple_bams_itd(
            filepaths,
            work_dir=work_dir,
            sample_id=sample_id,
            target_panel=target_panel,
            hotspots=hotspots,
            batch_size=batch_size,
            log=log,
        )
        job.context.add_metadata("itd_analysis", result)
        if result.get("files_processed", 0) == 0 and result.get("total_files", 0) > 0:
            job.context.add_error("itd_analysis", "No BAM files staged successfully")
        else:
            job.context.add_result(
                "itd_analysis",
                {
                    "success": True,
                    "sample_id": sample_id,
                    "events": result.get("events", 0),
                    "events_path": result.get("events_path"),
                    "genes": result.get("genes", []),
                },
            )
        return

    # Single-file path
    bam_path = job.context.filepath
    sample_id = _resolve_sample_id(job)
    job_panel = job.context.metadata.get("target_panel", target_panel)
    if job_panel != target_panel:
        log.warning(
            "Panel mismatch: metadata=%s handler=%s; using metadata",
            job_panel,
            target_panel,
        )
        target_panel = job_panel
        try:
            hotspots = _load_active_hotspots(
                target_panel, job.context.metadata, work_dir=work_dir
            )
        except Exception as exc:
            job.context.add_error("itd_analysis", str(exc))
            return

    log.info(
        "ITD analysis for %s (sample=%s, panel=%s, %d genes)",
        bam_path,
        sample_id,
        target_panel,
        len(hotspots),
    )

    try:
        stats, _should_accumulate = process_bam_itd_with_staging(
            bam_path,
            work_dir,
            sample_id,
            hotspots,
            batch_size=batch_size,
        )
        # Always rebuild reports after each file so the GUI stays current during
        # adaptive sampling; staging still appends dataset parts.
        accumulation = accumulate_itd_candidates(
            work_dir,
            sample_id,
            hotspots,
            force=True,
            batch_size=batch_size,
        )
        result = {
            "success": True,
            "sample_id": sample_id,
            "rows": stats.get("rows", 0),
            "pending": stats.get("pending"),
            "events": accumulation.get("events", 0),
            "events_path": accumulation.get("events_path"),
            "genes": sorted(hotspots),
        }
        job.context.add_metadata("itd_analysis", result)
        job.context.add_result("itd_analysis", result)
        log.info(
            "ITD complete for %s: %d events",
            sample_id,
            result.get("events", 0),
        )
    except Exception as exc:
        log.exception("ITD analysis failed: %s", exc)
        job.context.add_error("itd_analysis", str(exc))
