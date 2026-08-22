"""
Workflow hooks for integrating with the GUI message queue system.

This module provides hooks that send real-time updates to the GUI
without blocking workflow execution.
"""

import logging
import time
from typing import Any, Dict, List

from .gui_launcher import UpdateType, send_gui_update


def install_workflow_hooks(
    workflow_runner: Any, workflow_steps: List[str], monitored_directory: str
) -> None:
    """Install hooks into the workflow system for GUI updates.

    The hooks are intentionally non-invasive: we wrap public methods and start a
    lightweight polling thread that reads `manager.get_stats()` to drive the UI.
    """
    try:
        logging.info("Installing workflow hooks for GUI integration")

        # Prime the UI with an initial status
        try:
            send_gui_update(
                UpdateType.WORKFLOW_STATUS,
                {
                    "is_running": True,
                    "start_time": time.time(),
                    "workflow_steps": workflow_steps,
                    "monitored_directory": monitored_directory,
                },
                priority=10,
            )
        except Exception:
            pass

        _install_runner_hooks(workflow_runner, workflow_steps, monitored_directory)
        _install_manager_hooks(workflow_runner)

        logging.info("Workflow hooks installed for real-time GUI monitoring")
    except Exception as exc:
        logging.error(f"Failed to install workflow hooks: {exc}")
        logging.info("GUI will run without real-time updates")


def _infer_sample_id_from_path(filepath: str) -> str:
    """Best-effort sample id inference from a BAM path.

    Prefer the nearest ancestor directory whose name is not a BAM/FASTQ bucket
    (e.g., bam_fail, bam_pass, bam) and is not simply 'pass'/'fail'.
    Falls back to the immediate parent directory name.
    """
    try:
        from pathlib import Path

        p = Path(filepath)
        disallowed = {
            "bam_fail",
            "bam_pass",
            "bam",
            "fastq",
            "fastq_fail",
            "fastq_pass",
            "pass",
            "fail",
        }
        # Look up to 5 ancestors for a reasonable sample folder
        current = p.parent
        for _ in range(5):
            name = current.name
            if not name:
                break
            lname = name.lower()
            if lname not in disallowed and "bam" not in lname and "fastq" not in lname:
                return name
            current = current.parent
        return p.parent.name or "unknown"
    except Exception:
        return "unknown"


def _install_runner_hooks(
    workflow_runner: Any, workflow_steps: List[str], monitored_directory: str
) -> None:
    """Wrap `run_workflow` to emit start/stop UI events."""
    try:
        if not hasattr(workflow_runner, "run_workflow"):
            return

        original_run = workflow_runner.run_workflow

        def hooked_run_workflow(*args, **kwargs):
            logging.info("[GUI] run_workflow wrapper: start event will be sent")
            send_gui_update(
                UpdateType.WORKFLOW_STATUS,
                {
                    "is_running": True,
                    "start_time": time.time(),
                    "workflow_steps": workflow_steps,
                    "monitored_directory": monitored_directory,
                },
                priority=10,
            )
            try:
                return original_run(*args, **kwargs)
            finally:
                logging.info("[GUI] run_workflow wrapper: stop event will be sent")
                send_gui_update(
                    UpdateType.WORKFLOW_STATUS,
                    {"is_running": False, "stop_time": time.time()},
                    priority=10,
                )

        workflow_runner.run_workflow = hooked_run_workflow
    except Exception as exc:
        logging.debug(f"Failed to install runner hooks: {exc}")


def _install_manager_hooks(workflow_runner: Any) -> None:
    """Start a polling thread based on the runner's manager."""
    try:
        manager = getattr(workflow_runner, "manager", None) or getattr(
            workflow_runner, "workflow_manager", None
        )
        if manager is None:
            logging.debug("No workflow manager found on runner; skipping manager hooks")
            return
        logging.info("[GUI] Starting polling thread for GUI updates")
        _start_polling_updates(manager)
    except Exception as exc:
        logging.debug(f"Failed to install manager hooks: {exc}")


def _start_polling_updates(manager: Any, interval_seconds: float = 15.0) -> None:
    """Start a background thread that polls `manager.get_stats()` and sends GUI updates.

    Reduced default interval from 10s to 15s to reduce update frequency and prevent queue buildup.
    """
    import threading

    def poll() -> None:
        last_log_sent = 0.0
        while True:
            try:
                start_poll = time.time()
                stats = None
                try:
                    stats = manager.get_stats()
                except Exception as inner_exc:
                    logging.warning(f"[GUI] Poll get_stats() failed: {inner_exc}")
                    # continue to next iteration after sleep
                    raise
                finally:
                    elapsed = time.time() - start_poll
                    if elapsed > 5.0:
                        logging.warning(f"[GUI] Poll slow get_stats: {elapsed:.2f}s")
                now = time.time()
                logging.info(
                    f"[GUI] Poll: active={stats.get('active_jobs', 0)}, completed={stats.get('completed', 0)}, "
                    f"failed={stats.get('failed', 0)}, total={stats.get('total_actual_jobs', 0)}"
                )

                # Queue status (support legacy manager and new Ray-Coordinator stats)
                try:
                    # 1) Totals
                    submitted_by_queue = stats.get("submitted_by_queue", {}) or {}
                    totals_by_category = stats.get("totals_by_category", {}) or {}
                    queue_sizes = stats.get("queue_sizes", {}) or {}

                    def _total_for(q: str) -> int:
                        # Prefer legacy submitted_by_queue if present
                        if submitted_by_queue:
                            return int(submitted_by_queue.get(q, 0) or 0)
                        # Map new coordinator categories to the same keys
                        cat_map = {
                            "preprocessing": "preprocessing",
                            "mgmt": "mgmt",
                            "cnv": "cnv",
                            "target": "target",
                            "fusion": "fusion",
                            "classification": "classification",
                            "other": "other",
                            "slow": "other",
                        }
                        if totals_by_category:
                            key = cat_map.get(q, q)
                            return int(totals_by_category.get(key, 0) or 0)
                        # Fallback: instantaneous queue sizes if nothing else available
                        # Map 'other' to 'slow' in legacy queue_sizes
                        if q == "other":
                            return int(queue_sizes.get("slow", 0) or 0)
                        return int(queue_sizes.get(q, 0) or 0)

                    queue_update: Dict[str, Dict[str, int]] = {
                        "preprocessing": {
                            "running": 0,
                            "total": _total_for("preprocessing"),
                        },
                        "mgmt": {"running": 0, "total": _total_for("mgmt")},
                        "cnv": {"running": 0, "total": _total_for("cnv")},
                        "target": {"running": 0, "total": _total_for("target")},
                        "fusion": {"running": 0, "total": _total_for("fusion")},
                        "classification": {
                            "running": 0,
                            "total": _total_for("classification"),
                        },
                        "other": {"running": 0, "total": _total_for("other")},
                    }

                    # 2) Running
                    running_by_category = stats.get("running_by_category", {}) or {}
                    if running_by_category:
                        # Direct consumption from new coordinator
                        for k in [
                            "preprocessing",
                            "mgmt",
                            "cnv",
                            "target",
                            "fusion",
                            "classification",
                            "other",
                        ]:
                            queue_update[k]["running"] = int(
                                running_by_category.get(k, 0) or 0
                            )
                    else:
                        # Legacy path: derive from active_by_worker prefixes OR from active_by_queue lengths
                        active_by_worker = stats.get("active_by_worker", {}) or {}
                        if active_by_worker:
                            for worker_name, jobs in active_by_worker.items():
                                wn = str(worker_name)
                                n = len(jobs) if isinstance(jobs, list) else 0
                                if (
                                    wn.startswith("PreprocessingWorker")
                                    or wn == "preprocessing"
                                ):
                                    queue_update["preprocessing"]["running"] += n
                                elif wn.startswith("MGMTWorker") or wn == "mgmt":
                                    queue_update["mgmt"]["running"] += n
                                elif wn.startswith("CNVWorker") or wn == "cnv":
                                    queue_update["cnv"]["running"] += n
                                elif wn.startswith("TargetWorker") or wn == "target":
                                    queue_update["target"]["running"] += n
                                elif wn.startswith("FusionWorker") or wn == "fusion":
                                    queue_update["fusion"]["running"] += n
                                elif (
                                    wn.startswith("ClassificationWorker")
                                    or wn == "classification"
                                ):
                                    queue_update["classification"]["running"] += n
                                elif wn.startswith("SlowWorker") or wn in {
                                    "other",
                                    "slow",
                                }:
                                    queue_update["other"]["running"] += n

                    send_gui_update(UpdateType.QUEUE_UPDATE, queue_update, priority=4)
                except Exception:
                    pass

                # Active jobs table
                try:
                    job_rows: List[Dict[str, Any]] = []
                    for worker_name, jobs in stats.get("active_by_worker", {}).items():  # type: ignore
                        for job in jobs:
                            job_rows.append(
                                {
                                    "job_id": "",
                                    "job_type": job.get("job_type", ""),
                                    "filepath": job.get("filepath", ""),
                                    "worker_name": worker_name,
                                    "duration": int(job.get("duration", 0)),
                                    "progress": 0.0,
                                }
                            )
                    send_gui_update(
                        UpdateType.JOB_UPDATE, {"active_jobs": job_rows}, priority=5
                    )
                except Exception:
                    pass

                # Overall progress (include counts for the UI summary)
                try:
                    total_processed = stats.get("total_processed", 0)
                    total_actual_jobs = stats.get("total_actual_jobs", 0)
                    progress = (
                        min(0.999, float(total_processed) / float(total_actual_jobs))
                        if total_actual_jobs > 0
                        else 0.0
                    )
                    send_gui_update(
                        UpdateType.PROGRESS_UPDATE,
                        {
                            "progress": progress,
                            "completed": stats.get("completed", 0),
                            "failed": stats.get("failed", 0),
                            "total": total_actual_jobs,
                        },
                        priority=2,
                    )
                except Exception:
                    pass

                # Periodic log heartbeat (every ~20s) so users see activity in GUI
                if now - last_log_sent > 20:
                    try:
                        msg = (
                            f"Heartbeat: active={stats.get('active_jobs', 0)}, "
                            f"completed={stats.get('completed', 0)}, failed={stats.get('failed', 0)}"
                        )
                        send_gui_update(
                            UpdateType.LOG_MESSAGE,
                            {"level": "INFO", "message": msg},
                            priority=1,
                        )
                        last_log_sent = now
                    except Exception:
                        pass

                # Samples overview using manager-provided per-sample stats when available
                try:
                    if "samples" in stats and isinstance(stats["samples"], list):
                        send_gui_update(
                            UpdateType.SAMPLES_UPDATE,
                            {"samples": stats["samples"]},
                            priority=1,
                        )
                    else:
                        # Fallback to aggregation from active_by_worker (legacy)
                        samples_map: Dict[str, Dict[str, Any]] = {}
                        for _worker, jobs in stats.get("active_by_worker", {}).items():  # type: ignore
                            for job in jobs:
                                fp = job.get("filepath", "")
                                sid = job.get(
                                    "sample_id"
                                ) or _infer_sample_id_from_path(fp)
                                if sid.lower() in {
                                    "bam_fail",
                                    "bam_pass",
                                    "bam",
                                    "pass",
                                    "fail",
                                }:
                                    sid = _infer_sample_id_from_path(fp)
                                entry = samples_map.setdefault(
                                    sid,
                                    {
                                        "sample_id": sid,
                                        "active_jobs": 0,
                                        "total_jobs": 0,
                                        "completed_jobs": 0,
                                        "failed_jobs": 0,
                                        "job_types": set(),
                                        "last_seen": now,
                                    },
                                )
                                entry["active_jobs"] += 1
                                jt = job.get("job_type", "")
                                if jt:
                                    entry["job_types"].add(jt)
                                entry["last_seen"] = now
                        global_completed = stats.get("completed", 0)
                        global_failed = stats.get("failed", 0)
                        global_total = stats.get("total_actual_jobs", 0)
                        samples_payload = []
                        for entry in samples_map.values():
                            samples_payload.append(
                                {
                                    "sample_id": entry["sample_id"],
                                    "active_jobs": entry["active_jobs"],
                                    "total_jobs": global_total,
                                    "completed_jobs": global_completed,
                                    "failed_jobs": global_failed,
                                    "job_types": sorted(entry["job_types"]),
                                    "last_seen": entry["last_seen"],
                                }
                            )
                        send_gui_update(
                            UpdateType.SAMPLES_UPDATE,
                            {"samples": samples_payload},
                            priority=1,
                        )
                except Exception:
                    pass

            except Exception as exc:
                logging.error(f"[GUI] Polling error: {exc}")
            finally:
                time.sleep(interval_seconds)

    thread = threading.Thread(target=poll, daemon=True, name="robin-GUI-Polling")
    thread.start()
