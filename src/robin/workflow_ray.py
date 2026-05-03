"""
Ray Core implementation of the robin workflow engine
- Specialized per-queue actors (preprocessing, bed_conversion, mgmt, cnv, target, fusion, classification, slow)
- Central Coordinator actor for dedup (1 running + 1 pending per (sample_id, job_type)), triggers, stats
- Handlers registered for real robin job types, with resource hints per job
- tqdm-based live monitor similar to the original implementation

Usage (example):

    pip install "ray>=2.30" tqdm watchdog

    python ray_robin_core.py \
        --plan preprocessing:preprocessing \
        --plan bed_conversion:bed_conversion \
        --plan mgmt:mgmt \
        --plan sturgeon:sturgeon \
        --paths /data/incoming

Replace/extend the plan and paths as needed. If you already have a Watchdog-based
file watcher, call submit_jobs() on the Coordinator in your event handlers instead.
"""

from __future__ import annotations

import logging

# Suppress pkg_resources deprecation warnings from sorted_nearest
import warnings
warnings.filterwarnings(
    "ignore", message="pkg_resources is deprecated", category=UserWarning
)
# Suppress matplotlib tight_layout warnings
warnings.filterwarnings(
    "ignore", message="The figure layout has changed to tight", category=UserWarning
)

import os
import threading
import time
import tempfile
import pickle
import uuid
import asyncio
import argparse
from collections import deque
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple, Callable
import inspect
from pathlib import Path
import itertools
from contextlib import nullcontext

import ray
from tqdm import tqdm
from watchdog.events import FileSystemEventHandler
from watchdog.observers import Observer

try:
    from rich.progress import (
        BarColumn,
        Progress,
        SpinnerColumn,
        TaskProgressColumn,
        TextColumn,
        TimeElapsedColumn,
        TimeRemainingColumn,
    )
    from rich.console import Console

    _RICH_AVAILABLE = True
except Exception:
    _RICH_AVAILABLE = False

_RICH_CONSOLE = Console() if _RICH_AVAILABLE else None


def _print_styled(message: str, level: str = "info") -> None:
    """Print with optional rich styling (no emojis)."""
    if not _RICH_AVAILABLE or _RICH_CONSOLE is None:
        print(message)
        return

    styles = {
        "info": "cyan",
        "success": "green",
        "warn": "yellow",
        "error": "red",
        "header": "bold magenta",
        "muted": "dim",
    }
    style = styles.get(level, "white")
    _RICH_CONSOLE.print(message, style=style)

# Import memory management
try:
    from robin.memory_manager import MemoryManager
except ImportError:
    MemoryManager = None

# Disable memory manager for testing via env flag
DISABLE_MEMORY_MANAGER = True #os.getenv("ROBIN_DISABLE_MEMORY_MANAGER", "0") == "1"
_HANDLER_IMPORT_ERRORS: Dict[str, str] = {}

# Optional GUI hook integration
try:
    from robin.gui_launcher import (
        send_gui_update as _gui_send_update,
        UpdateType as _GUIUpdateType,
        launch_gui as _gui_launch,
    )
except Exception as e:
    raise Exception(f"GUI not available: {e}")

    def _gui_send_update(*args, **kwargs):
        return None

    _GUIUpdateType = None

    def _gui_launch(*args, **kwargs):
        return None


# Import handlers
try:
    from robin.analysis.bam_preprocessor import (
        bam_preprocessing_handler as _preprocessing_handler,
    )
except Exception:
    _preprocessing_handler = None

try:
    from robin.analysis.bed_conversion import (
        bed_conversion_handler as _bed_conversion_handler,
    )
except Exception:
    _bed_conversion_handler = None

try:
    from robin.analysis.mgmt_analysis import mgmt_handler as _mgmt_handler
except Exception:
    _mgmt_handler = None

try:
    from robin.analysis.cnv_analysis import (
        cnv_handler as _cnv_handler,
        clear_sample_cache,
        get_sample_cache_stats,
        cleanup_sample_cache_on_completion
    )
except Exception as exc:
    _cnv_handler = None
    clear_sample_cache = None
    get_sample_cache_stats = None
    cleanup_sample_cache_on_completion = None
    _HANDLER_IMPORT_ERRORS["cnv"] = str(exc)

try:
    from robin.analysis.target_analysis import (
        target_handler as _target_handler,
        igv_bam_handler as _igv_bam_handler,
        snp_analysis_handler as _snp_analysis_handler,
        target_bam_finalize_handler as _target_bam_finalize_handler,
    )
except Exception:
    _target_handler = None
    _igv_bam_handler = None
    _snp_analysis_handler = None
    _target_bam_finalize_handler = None

try:
    from robin.analysis.fusion_analysis import fusion_handler as _fusion_handler
except Exception as exc:
    _fusion_handler = None
    _HANDLER_IMPORT_ERRORS["fusion"] = str(exc)

try:
    from robin.analysis.sturgeon_analysis import (
        sturgeon_handler as _sturgeon_handler,
    )
except Exception:
    _sturgeon_handler = None

try:
    from robin.analysis.nanodx_analysis import nanodx_handler as _nanodx_handler
except Exception:
    _nanodx_handler = None

try:
    from robin.analysis.nanodx_analysis import (
        pannanodx_handler as _pannanodx_handler,
    )
except Exception:
    _pannanodx_handler = None

try:
    from robin.analysis.random_forest_analysis import (
        random_forest_handler as _random_forest_handler,
    )
except Exception:
    _random_forest_handler = None

# Optional logging helper
try:
    from robin.logging_config import (
        get_job_logger as _get_job_logger,
        configure_logging as _configure_logging,
    )
except Exception:

    def _get_job_logger(job_id: str, job_type: str, filepath: str):
        class _L:
            def info(self, *a, **k):
                pass

            def debug(self, *a, **k):
                pass

            def error(self, *a, **k):
                pass

        return _L()

    def _configure_logging(
        global_level: str = "INFO", job_levels: Optional[Dict[str, str]] = None
    ):
        return None


GLOBAL_LOG_LEVEL: str = "INFO"


# ---------- Batch Configuration ----------
BATCH_CONFIG: Dict[str, Dict[str, Any]] = {
    # Preprocessing should NOT be batched - each file needs individual sample ID extraction
    "preprocessing": {"max_batch_size": 1, "timeout_seconds": 0},  # Force no batching
    "bed_conversion": {"max_batch_size": 20, "timeout_seconds": 2},
    "mgmt": {"max_batch_size": 20, "timeout_seconds": 2},
    "cnv": {"max_batch_size": 20, "timeout_seconds": 2},
    "target": {"max_batch_size": 20, "timeout_seconds": 2},
    "fusion": {"max_batch_size": 20, "timeout_seconds": 2},
    "sturgeon": {"max_batch_size": 20, "timeout_seconds": 2},
    "nanodx": {"max_batch_size": 20, "timeout_seconds": 2},
    "pannanodx": {"max_batch_size": 20, "timeout_seconds": 2},
    "random_forest": {"max_batch_size": 20, "timeout_seconds": 2},
    "igv_bam": {"max_batch_size": 20, "timeout_seconds": 2},
    "snp_analysis": {"max_batch_size": 20, "timeout_seconds": 2},
}

# ---------- Shared DTOs ----------
@dataclass
class WorkflowContext:
    filepath: str
    metadata: Dict[str, Any] = field(default_factory=dict)
    results: Dict[str, Any] = field(default_factory=dict)
    history: List[str] = field(default_factory=list)
    errors: List[Dict[str, Any]] = field(default_factory=list)
    
    # NEW: Simple batch metadata
    batch_id: Optional[str] = None
    batch_index: Optional[int] = None

    def add_metadata(self, key: str, value: Any) -> None:
        self.metadata[key] = value

    def add_result(self, job_type: str, result: Any) -> None:
        self.results[job_type] = result
        self.history.append(job_type)

    def add_error(self, job_type: str, error: str) -> None:
        self.errors.append(
            {"job_type": job_type, "error": error, "timestamp": time.time()}
        )

    def get_sample_id(self) -> str:
        # First try to get from bam_metadata (for backward compatibility)
        bam_md = self.metadata.get("bam_metadata", {})
        sample_id = bam_md.get("sample_id", "unknown")
        
        # If not found or "unknown", try to get from preprocessing results
        if sample_id == "unknown":
            preprocessing_result = self.results.get("preprocessing", {})
            sample_id = preprocessing_result.get("sample_id", "unknown")
        
        return sample_id
    
    def set_batch_info(self, batch_id: str, batch_index: int) -> None:
        """Set batch information for this context"""
        self.batch_id = batch_id
        self.batch_index = batch_index
        self.metadata["batch_id"] = batch_id
        self.metadata["batch_index"] = batch_index
    
    def cleanup_intermediate_data(self) -> None:
        """
        Clean up intermediate data from context to reduce memory usage.
        This removes large data structures that are no longer needed after processing.
        """
        # Remove large lists from results that are stored on disk
        for job_type in list(self.results.keys()):
            result = self.results.get(job_type)
            if isinstance(result, dict):
                # Remove large lists from preprocessing results
                if job_type == "preprocessing" and "supplementary_read_ids" in result:
                    result["supplementary_read_ids"] = []  # Clear large list
                    result["supplementary_read_ids_cleared"] = True
        
        # Keep history but limit errors list to last 10 errors
        if len(self.errors) > 10:
            self.errors = self.errors[-10:]


@dataclass
class Job:
    job_id: int
    job_type: str
    origin: str
    workflow: List[str]
    step: int
    context: WorkflowContext

    def next_job(self) -> Optional["Job"]:
        if self.step + 1 >= len(self.workflow):
            return None
        next_step = self.workflow[self.step + 1]
        if ":" not in next_step:
            raise ValueError(f"Invalid workflow step: {next_step}")
        queue_type, next_type = next_step.split(":", 1)
        return Job(
            self.job_id,
            next_type,
            queue_type,
            self.workflow,
            self.step + 1,
            self.context,
        )


@dataclass
class BatchedJob:
    job_id: int
    job_type: str
    origin: str
    workflow: List[str]
    step: int
    contexts: List[WorkflowContext]  # Multiple contexts for batched processing
    batch_id: str
    sample_id: str
    
    def get_sample_id(self) -> str:
        return self.sample_id
    
    def get_file_count(self) -> int:
        return len(self.contexts)
    
    def get_filepaths(self) -> List[str]:
        return [ctx.filepath for ctx in self.contexts]


# ---------- Queue mapping & triggers ----------
QUEUE_TO_TYPES: Dict[str, Set[str]] = {
    "preprocessing": {"preprocessing"},
    "bed_conversion": {"bed_conversion"},
    "mgmt": {"mgmt"},
    "cnv": {"cnv"},
    "target": {"target"},
    "fusion": {"fusion"},
    "classification": {"sturgeon", "nanodx", "pannanodx"},
    "slow": {"random_forest", "igv_bam", "snp_analysis", "target_bam_finalize"},
}

TRIGGERS: Dict[str, List[str]] = {
    # preprocessing -> analyses
    "preprocessing": ["bed_conversion", "mgmt", "cnv", "target", "fusion"],
    # bed_conversion -> classifiers
    "bed_conversion": ["sturgeon", "nanodx", "pannanodx", "random_forest"],
    # Build IGV-ready BAM after target analysis
    "target": ["igv_bam"],
}

# Per-sample de-duplication to avoid output races. Ensure only one job of these types
# runs concurrently per sample (max 1 running + 1 pending).
# Only deduplicate selected per-sample jobs; analysis/bed_conversion should queue, not skip.
# Include SNP analysis to avoid duplicate manual submissions for the same sample.
DEDUP_TYPES: Set[str] = {
    "sturgeon",
    "nanodx",
    "pannanodx",
    "random_forest",
    "snp_analysis",
}

# Per-sample serialization by type to avoid races on shared per-sample outputs
# Disabled: treat CNV like other analysis jobs (mgmt/target/fusion)
SERIALIZE_BY_TYPE_PER_SAMPLE: Set[str] = set()

# Classification job types (single global pipeline per type)
CLASSIFICATION_TYPES: Set[str] = {"sturgeon", "nanodx", "pannanodx", "random_forest"}

# Job types that must be serialized per sample (no overlap across these types)
## Removed per-sample cross-type serialization to allow one job per type globally

# Ray memory is in bytes. Default 1 GiB per task; fusion 4 GiB; snp_analysis 8 GiB (Clair3/variant calling).
_GB = 1024 * 1024 * 1024
_DEFAULT_MEMORY = 1 * _GB
_FUSION_MEMORY = 4 * _GB
_SNP_ANALYSIS_MEMORY = 8 * _GB
RESOURCE_HINTS: Dict[str, Dict[str, Any]] = {
    # Tune these to your cluster
    "preprocessing": {"num_cpus": 1, "memory": _DEFAULT_MEMORY},
    "bed_conversion": {"num_cpus": 1, "memory": _DEFAULT_MEMORY},
    "mgmt": {"num_cpus": 1, "memory": _DEFAULT_MEMORY},
    "cnv": {"num_cpus": 1, "memory": _DEFAULT_MEMORY},  # CNV gets dedicated CPU but only 1 thread
    "target": {"num_cpus": 1, "memory": _DEFAULT_MEMORY},
    "fusion": {"num_cpus": 1, "memory": _FUSION_MEMORY},  # Fusion needs 4 GiB (BAM scan, breakpoint aggregation)
    # Classifiers do not require GPU by default (CPU-only)
    "sturgeon": {"num_cpus": 1, "memory": _DEFAULT_MEMORY},
    "nanodx": {"num_cpus": 1, "memory": _DEFAULT_MEMORY},
    "pannanodx": {"num_cpus": 1, "memory": _DEFAULT_MEMORY},
    "random_forest": {"num_cpus": 1, "memory": _DEFAULT_MEMORY},
    "igv_bam": {"num_cpus": 1, "memory": _DEFAULT_MEMORY},
    "snp_analysis": {"num_cpus": 1, "memory": _SNP_ANALYSIS_MEMORY},  # Clair3/variant calling needs 8 GiB
    "target_bam_finalize": {"num_cpus": 1, "memory": _DEFAULT_MEMORY},
}

# ---------- Sample Job Batcher ----------
class SampleJobBatcher:
    def __init__(self):
        # sample_id -> job_type -> List[Job]
        self.pending_jobs: Dict[str, Dict[str, List[Job]]] = {}
        # sample_id -> job_type -> timestamp
        self.last_job_time: Dict[str, Dict[str, float]] = {}
        # Lock created lazily to avoid serialization issues in Ray actors
        self._lock: Optional[asyncio.Lock] = None
    
    def _get_lock(self) -> Optional[asyncio.Lock]:
        """Get or create the asyncio lock lazily."""
        if self._lock is None and asyncio:
            self._lock = asyncio.Lock()
        return self._lock
    
    async def add_job(self, job: Job) -> List[BatchedJob]:
        """Add a job and return any completed batches"""
        sample_id = job.context.get_sample_id()
        job_type = job.job_type
        
        lock = self._get_lock()
        async with lock if lock else nullcontext():
            # Initialize if needed
            if sample_id not in self.pending_jobs:
                self.pending_jobs[sample_id] = {}
                self.last_job_time[sample_id] = {}
            
            if job_type not in self.pending_jobs[sample_id]:
                self.pending_jobs[sample_id][job_type] = []
                self.last_job_time[sample_id][job_type] = time.time()
            
            # Add job to pending list
            self.pending_jobs[sample_id][job_type].append(job)
            self.last_job_time[sample_id][job_type] = time.time()
            
            # Check for completed batches
            batches = self._check_and_create_batches(sample_id, job_type)
            return batches
    
    def _check_and_create_batches(self, sample_id: str, job_type: str) -> List[BatchedJob]:
        """Check if we should create batches for a sample/job_type combination.
        Jobs with force_individual_batch (e.g. large BAMs) are emitted as single-file batches.
        """
        config = BATCH_CONFIG.get(job_type, {"max_batch_size": 20, "timeout_seconds": 10})
        max_batch_size = config["max_batch_size"]
        
        pending = self.pending_jobs[sample_id][job_type]
        batches = []
        
        # Emit single-file batches for jobs marked force_individual_batch (e.g. large BAMs)
        individual = [j for j in pending if (j.context.metadata or {}).get("force_individual_batch")]
        rest = [j for j in pending if not (j.context.metadata or {}).get("force_individual_batch")]
        for job in individual:
            batches.append(self._create_batched_job([job], sample_id, job_type))
        
        # Create batches of max_batch_size from the rest
        while len(rest) >= max_batch_size:
            batch_jobs = rest[:max_batch_size]
            rest = rest[max_batch_size:]
            batched_job = self._create_batched_job(batch_jobs, sample_id, job_type)
            batches.append(batched_job)
        
        # Update pending list (only unbatched jobs remain)
        self.pending_jobs[sample_id][job_type] = rest
        
        return batches
    
    async def check_timeouts(self) -> List[BatchedJob]:
        """Check for timed-out batches and return them"""
        current_time = time.time()
        timed_out_batches = []
        
        lock = self._get_lock()
        async with lock if lock else nullcontext():
            for sample_id in list(self.pending_jobs.keys()):
                for job_type in list(self.pending_jobs[sample_id].keys()):
                    jobs = self.pending_jobs[sample_id][job_type]
                    if not jobs:
                        continue
                    
                    config = BATCH_CONFIG.get(job_type, {"timeout_seconds": 10})
                    timeout_seconds = config["timeout_seconds"]
                    last_time = self.last_job_time[sample_id][job_type]
                    
                    if (current_time - last_time) >= timeout_seconds:
                        # Create batch with remaining jobs
                        batch = self._create_batched_job(jobs, sample_id, job_type)
                        timed_out_batches.append(batch)
                        
                        # Clear the pending jobs
                        self.pending_jobs[sample_id][job_type] = []
                        del self.last_job_time[sample_id][job_type]
            
            # Clean up empty entries
            self._cleanup_empty_entries()
        
        return timed_out_batches
    
    def _create_batched_job(self, jobs: List[Job], sample_id: str, job_type: str) -> BatchedJob:
        """Create a batched job from a list of individual jobs"""
        if not jobs:
            raise ValueError("Cannot create batched job from empty job list")
        
        # Validate all jobs have same sample_id and job_type
        for job in jobs:
            if job.context.get_sample_id() != sample_id:
                raise ValueError(f"Mixed sample IDs in batch: {sample_id} vs {job.context.get_sample_id()}")
            if job.job_type != job_type:
                raise ValueError(f"Mixed job types in batch: {job_type} vs {job.job_type}")
        
        # Use the first job as template
        template_job = jobs[0]
        batch_id = f"{sample_id}_{job_type}_{int(time.time() * 1000)}"
        
        # Extract contexts from all jobs
        contexts = [job.context for job in jobs]
        
        return BatchedJob(
            job_id=next(_job_id_counter),
            job_type=job_type,
            origin=template_job.origin,
            workflow=template_job.workflow,
            step=template_job.step,
            contexts=contexts,
            batch_id=batch_id,
            sample_id=sample_id
        )
    
    def _cleanup_empty_entries(self):
        """Remove empty entries from pending jobs and timestamps"""
        # Remove empty job type entries
        for sample_id in list(self.pending_jobs.keys()):
            for job_type in list(self.pending_jobs[sample_id].keys()):
                if not self.pending_jobs[sample_id][job_type]:
                    del self.pending_jobs[sample_id][job_type]
                    if job_type in self.last_job_time[sample_id]:
                        del self.last_job_time[sample_id][job_type]
            
            # Remove empty sample entries
            if not self.pending_jobs[sample_id]:
                del self.pending_jobs[sample_id]
                if sample_id in self.last_job_time:
                    del self.last_job_time[sample_id]


# ---------- Utilities ----------
_job_id_counter = itertools.count(1000)


def job_queue_of(job_type: str) -> str:
    for q, types in QUEUE_TO_TYPES.items():
        if job_type in types:
            return q
    return "slow"


def _get_ray_runtime_ids() -> Tuple[Optional[str], Optional[str], Optional[str]]:
    try:
        ctx = ray.get_runtime_context()
        task_id = str(ctx.get_task_id())
        job_id = str(ctx.get_job_id())
        node_id = None
        try:
            node_id = str(ctx.get_node_id())
        except Exception:
            node_id = None
        return task_id, job_id, node_id
    except Exception:
        return None, None, None


def _get_sample_id(job: Job) -> Optional[str]:
    try:
        if hasattr(job.context, "get_sample_id"):
            sample_id = job.context.get_sample_id()
            if sample_id:
                return sample_id
    except Exception:
        pass
    try:
        return job.context.metadata.get("sample_id")
    except Exception:
        return None


# ---------- Real handler wrappers as Ray tasks ----------
# Each wrapper returns an updated WorkflowContext.

NEEDS_WORK_DIR: Set[str] = {
    "bed_conversion",
    "mgmt",
    "cnv",
    "target",
    "fusion",
    "sturgeon",
    "nanodx",
    "pannanodx",
    "random_forest",
}


def _notify_coordinator_handler_started(job_id: int) -> None:
    """Tell Coordinator the Ray worker has begun running the handler (not queued in Pool)."""
    try:
        ray.get_actor("robin_coordinator").mark_handler_started.remote(job_id)
    except Exception:
        pass


def _wrap_real_handler(
    py_handler: Optional[Callable[[Job], None]], job_type: str
) -> Callable[[Job], WorkflowContext]:
    def _impl(job: Job) -> WorkflowContext:
        # Initialize memory manager for this handler execution
        memory_manager = None
        if MemoryManager is not None and not DISABLE_MEMORY_MANAGER:
            try:
                # Configure memory management based on job type
                gc_every = 10 if job_type in {"mgmt", "cnv", "target", "fusion"} else 25
                rss_trigger = 1024 if job_type in {"mgmt", "cnv", "target", "fusion"} else 1024
                memory_manager = MemoryManager(
                    gc_every=gc_every,
                    rss_trigger_mb=rss_trigger,
                    enable_malloc_trim=True
                )
            except Exception:
                # Fallback if memory manager fails to initialize
                pass
        
        # Ensure per-process logging honors global level
        try:
            _configure_logging(global_level=GLOBAL_LOG_LEVEL)
        except Exception:
            pass
        logger = _get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
        task_id, ray_job_id, node_id = _get_ray_runtime_ids()
        sample_id = _get_sample_id(job)
        filepath = getattr(job.context, "filepath", None)
        logger.info(
            f"ray_task_start job_type={job_type} "
            f"job_id={getattr(job, 'job_id', None)} "
            f"ray_job_id={ray_job_id} task_id={task_id} node_id={node_id} "
            f"pid={os.getpid()} sample_id={sample_id} filepath={filepath}"
        )
        # Progress UI "duration" prefers handler wall time (excludes Pool/Ray queue after dispatch).
        try:
            jid = getattr(job, "job_id", None)
            if jid is not None:
                _notify_coordinator_handler_started(int(jid))
        except Exception:
            pass
        try:
            if py_handler is None:
                import_reason = _HANDLER_IMPORT_ERRORS.get(job_type)
                if import_reason:
                    raise RuntimeError(
                        f"No implementation available for job_type='{job_type}'. "
                        f"Handler import failed: {import_reason}"
                    )
                raise RuntimeError(
                    f"No implementation available for job_type='{job_type}'"
                )
            # Call user's handler; if they mutate job.context in-place, we just return it
            work_dir = job.context.metadata.get("work_dir")
            # Ensure output directories exist when a work_dir is provided
            if work_dir:
                try:
                    os.makedirs(work_dir, exist_ok=True)
                except Exception:
                    pass
                try:
                    sid = (
                        job.context.get_sample_id()
                        if hasattr(job.context, "get_sample_id")
                        else None
                    )
                    if sid and sid != "unknown":
                        os.makedirs(os.path.join(work_dir, sid), exist_ok=True)
                except Exception:
                    pass
            try:
                sig = inspect.signature(py_handler)
                # Check if handler accepts reference parameter
                accepts_reference = "reference" in sig.parameters
                accepts_target_panel = "target_panel" in sig.parameters

                # Get reference from job metadata if available
                reference = job.context.metadata.get("reference")
                target_panel = job.context.metadata.get("target_panel")
                if not target_panel:
                    raise ValueError(f"No target_panel found in job metadata for {job_type} handler")

                # Call handler with appropriate parameters
                if (
                    "work_dir" in sig.parameters
                    and work_dir
                    and (job_type in NEEDS_WORK_DIR)
                ):
                    if accepts_reference and reference and accepts_target_panel and job_type in ["fusion", "target", "cnv"]:
                        py_handler(job, work_dir=work_dir, reference=reference, target_panel=target_panel)
                    elif accepts_reference and reference and job_type in ["mgmt", "target"]:
                        py_handler(job, work_dir=work_dir, reference=reference)
                    elif accepts_target_panel and job_type in ["fusion", "target", "cnv"]:
                        py_handler(job, work_dir=work_dir, target_panel=target_panel)
                    else:
                        py_handler(job, work_dir=work_dir)
                elif accepts_reference and reference and accepts_target_panel and job_type in ["fusion", "target", "cnv"]:
                    py_handler(job, reference=reference, target_panel=target_panel)
                elif accepts_reference and reference and job_type in ["mgmt", "target"]:
                    py_handler(job, reference=reference)
                elif accepts_target_panel and job_type in ["fusion", "target", "cnv"]:
                    py_handler(job, target_panel=target_panel)
                else:
                    py_handler(job)
            except Exception:
                # Fallback to legacy call
                py_handler(job)
            # Log completion or skip based on handler-set context
            skipped = False
            try:
                skipped = bool(job.context.metadata.get("skip_downstream", False))
            except Exception:
                skipped = False
            if skipped and job_type == "preprocessing":
                logger.warning(
                    f"{job_type} skipped (deliberate): {os.path.basename(job.context.filepath)}"
                )
                try:
                    # Inform GUI/CLI via GUI hook if available
                    _gui_send_update(
                        _GUIUpdateType.LOG_MESSAGE,
                        {
                            "message": f"Skipped (too many reads): {os.path.basename(job.context.filepath)}",
                            "level": "WARNING",
                        },
                        priority=0,
                    )
                except Exception:
                    pass
            else:
                logger.info(f"{job_type} completed")
        except Exception as e:
            # Suppress expected downstream failures for fail-only BAM submissions.
            # Requirement: if *any* BAM in the submission contains 'pass', errors must remain visible.
            suppress = False
            try:
                if job_type in (CLASSIFICATION_TYPES | {"bed_conversion"}):
                    if bool(job.context.metadata.get("fail_only_bam_submission", False)):
                        fp = getattr(job.context, "filepath", "") or ""
                        if _bam_filename_kind(fp) == "fail":
                            suppress = True
            except Exception:
                suppress = False

            if suppress:
                # Expected: downstream tools may not support fail-only read BAMs.
                # We still attempted analysis; we just don't treat this as a failing job.
                try:
                    job.context.add_result(
                        job_type,
                        {
                            "status": "expected_failure",
                            "reason": "fail_only_bam_submission",
                            "message": str(e),
                        },
                    )
                except Exception:
                    pass
                logger.warning(
                    f"{job_type} expected failure for fail-only BAM submission: {e}"
                )
            else:
                job.context.add_error(job_type, str(e))
                logger.error(f"{job_type} failed: {e}")
                raise  # Re-raise so task fails and coordinator gets failure_kind="exception"
        finally:
            # Mark result if user handler didn't add one.
            # Keep a dict payload so downstream code can safely call .get(...).
            if job_type not in job.context.results:
                had_error = any(
                    (err.get("job_type") == job_type) for err in (job.context.errors or [])
                )
                fallback_result = {"status": "failed"} if had_error else {"status": "success"}
                job.context.add_result(job_type, fallback_result)
            
            # Trigger memory cleanup after handler execution
            # Let memory manager decide if cleanup is needed based on its heuristics
            if memory_manager is not None:
                try:
                    cleanup_stats = memory_manager.check_and_cleanup(is_idle=True)
                    # Log memory cleanup if it was triggered
                    if cleanup_stats.get('gc_triggered', False):
                        logger.debug(f"Memory cleanup: {cleanup_stats}")
                except Exception:
                    # Silently continue if memory management fails
                    pass
        
        return job.context

    return _impl


def enhanced_handler(job: BatchedJob) -> None:
    """Enhanced handler that processes batched jobs sequentially"""
    
    sample_id = job.get_sample_id()
    job_type = job.job_type
    batch_size = job.get_file_count()
    
    # Process each context sequentially
    for i, context in enumerate(job.contexts):
        context.set_batch_info(job.batch_id, i)
        
        try:
            # Process individual file within batch
            process_single_file_in_batch(context, job_type, i, batch_size)
            context.add_result(job_type, {"status": "success"})
            
        except Exception as e:
            # Fail entire batch if any file fails
            error_msg = f"File {i+1}/{batch_size} failed: {str(e)}"
            for ctx in job.contexts:
                ctx.add_error(job_type, error_msg)
            raise  # Re-raise to fail the entire batch
    
    # Mark batch as completed
    for context in job.contexts:
        context.add_result(job_type, {"status": "batch_completed"})

def process_single_file_in_batch(context: WorkflowContext, job_type: str, index: int, total: int) -> None:
    """Process a single file within a batch - to be implemented per job type"""
    # This would call the existing single-file processing logic
    pass


# Build remote versions with no fixed options; Pool will apply .options(**RESOURCE_HINTS[job_type])


# Real wrappers
preprocessing_handler_remote = ray.remote(
    _wrap_real_handler(_preprocessing_handler, "preprocessing")
)
bed_conversion_handler_remote = ray.remote(
    _wrap_real_handler(_bed_conversion_handler, "bed_conversion")
)
mgmt_handler_remote = ray.remote(_wrap_real_handler(_mgmt_handler, "mgmt"))
cnv_handler_remote = ray.remote(_wrap_real_handler(_cnv_handler, "cnv"))
target_handler_remote = ray.remote(_wrap_real_handler(_target_handler, "target"))
fusion_handler_remote = ray.remote(_wrap_real_handler(_fusion_handler, "fusion"))
sturgeon_handler_remote = ray.remote(_wrap_real_handler(_sturgeon_handler, "sturgeon"))
nanodx_handler_remote = ray.remote(_wrap_real_handler(_nanodx_handler, "nanodx"))
pannanodx_handler_remote = ray.remote(
    _wrap_real_handler(_pannanodx_handler, "pannanodx")
)
random_forest_handler_remote = ray.remote(
    _wrap_real_handler(_random_forest_handler, "random_forest")
)


# ---------- TypeProcessor & Coordinator Actors ----------
@ray.remote
class TypeProcessor:
    def __init__(self, job_type: str, remote_func, resource_options: Dict[str, Any]):
        self.job_type = job_type
        self.remote_func = remote_func
        self.resource_options = resource_options or {}
        
        # Initialize memory manager for long-running actor
        self.memory_manager = None
        if MemoryManager is not None and not DISABLE_MEMORY_MANAGER:
            try:
                # Configure memory management based on job type
                gc_every = 25 if job_type in {"mgmt", "cnv", "target", "fusion"} else 50
                rss_trigger = 1024 if job_type in {"mgmt", "cnv", "target", "fusion"} else 2048
                restart_every = 5000 if job_type in {"mgmt", "cnv", "target", "fusion"} else 10000
                restart_rss_trigger = 2048 if job_type in {"mgmt", "cnv", "target", "fusion"} else 4096
                self.memory_manager = MemoryManager(
                    gc_every=gc_every,
                    rss_trigger_mb=rss_trigger,
                    enable_malloc_trim=True,
                    restart_every=restart_every,
                    restart_rss_trigger_mb=restart_rss_trigger
                )
            except Exception as e:
                # Fallback if memory manager fails to initialize
                pass

    def update_handler(self, remote_func, resource_options: Dict[str, Any]):
        self.remote_func = remote_func
        self.resource_options = resource_options or {}
    
    def should_restart(self) -> bool:
        """Check if actor should restart due to memory management."""
        if self.memory_manager is not None:
            return self.memory_manager.is_restart_requested()
        return False
    
    def get_restart_reason(self) -> Optional[str]:
        """Get the reason for restart request."""
        if self.memory_manager is not None:
            return self.memory_manager.get_restart_reason()
        return None

    async def process(self, job: Job):
        # Submit the real handler task and await the result to avoid nested ObjectRefs
        ref = self.remote_func.options(**(self.resource_options or {})).remote(job)
        try:
            ctx = await ref
        except Exception:
            # Propagate error to coordinator; it will handle ok=False
            raise
        
        # Trigger memory cleanup after processing (actor is idle after job completion)
        if self.memory_manager is not None:
            try:
                cleanup_stats = self.memory_manager.check_and_cleanup(is_idle=True)
                # Log memory cleanup if it was triggered
                if cleanup_stats.get('gc_triggered', False):
                    pass  # Memory cleanup performed
                
                # Actor restart mechanism disabled - no restart checks
            except Exception as e:
                # Silently continue if memory management fails
                pass
        
        return ctx


def _job_timeout_seconds() -> int:
    """Per-job timeout in seconds (0 = disabled). From env ROBIN_JOB_TIMEOUT_SECONDS."""
    try:
        val = os.environ.get("ROBIN_JOB_TIMEOUT_SECONDS", "0").strip()
        return max(0, int(val))
    except (ValueError, TypeError):
        return 0


@ray.remote
class Coordinator:
    def __init__(
        self,
        target_panel: str,
        analysis_workers: int = 1,
        preset: Optional[str] = None,
        reference: Optional[str] = None,
        enable_batching: bool = True,
        job_timeout_seconds: int = 0,
    ):
        # dedup maps
        self.pending: Dict[Tuple[str, str], int] = {}
        self.running: Dict[Tuple[str, str], int] = {}
        # stats - use counters instead of lists to avoid unbounded growth
        self.completed_count: int = 0
        self.failed_count: int = 0
        self.completed_by_type: Dict[str, int] = {}
        self.failed_by_type: Dict[str, int] = {}
        self.submitted_by_type: Dict[str, int] = {}
        self.total_enqueued = 0
        self.total_skipped = 0
        self.active: Dict[int, Dict[str, Any]] = {}

        self.analysis_workers = analysis_workers
        self.target_panel = target_panel

        # Per-sample per-type serialization tracking
        self.running_by_type_sample: Dict[Tuple[str, str], int] = {}
        self.pending_by_type_sample: Dict[Tuple[str, str], int] = {}
        self.waiting_by_type_sample: Dict[Tuple[str, str], List[Job]] = {}

        # Per-sample one-shot scheduling for classifiers/slow
        self.scheduled_by_sample: Dict[str, Set[str]] = {}

        # Global one-in-flight (plus at most one waiting) per classification type
        self.classif_pending_by_type: Dict[str, int] = {}
        self.classif_waiting_by_type: Dict[str, Optional[Job]] = {}

        # create one TypeProcessor actor per job type
        self.processors: Dict[str, Any] = {}
        # inflight tracking: ObjectRef -> job
        self._inflight: Dict[Any, Job] = {}
        # backpressure: limit outstanding submissions per job type to avoid
        # massive pending tasks and worker explosion
        self.max_inflight_per_type: int = 64
        self.inflight_by_type: Dict[str, int] = {}
        self.waiting_by_type_global: Dict[str, List[Job]] = {}
        # backpressure: global and per-queue caps to prevent queue explosions
        queue_count = max(1, len(QUEUE_TO_TYPES))
        self.max_total_inflight: int = self.max_inflight_per_type * queue_count
        # Backpressure cap for total waiting jobs.
        # Test change: temporarily allow 2x waiting capacity to see if we
        # simply stall because the queue limit is reached.
        self.max_total_waiting: int = self.max_total_inflight * 40
        self.max_inflight_per_queue: Dict[str, int] = {
            q: self.max_inflight_per_type for q in QUEUE_TO_TYPES
        }
        self.max_waiting_per_queue: Dict[str, int] = {
            q: self.max_inflight_per_type * 200 for q in QUEUE_TO_TYPES
        }
        # Optional debug logging when backpressure caps are hit.
        # Enable with:
        #   ROBIN_DEBUG_BACKPRESSURE=1
        #   ROBIN_DEBUG_BACKPRESSURE_INTERVAL_S=10
        self._debug_backpressure: bool = os.getenv("ROBIN_DEBUG_BACKPRESSURE", "0") == "1"
        try:
            self._debug_backpressure_interval_s: float = float(
                os.getenv("ROBIN_DEBUG_BACKPRESSURE_INTERVAL_S", "10")
            )
        except Exception:
            self._debug_backpressure_interval_s = 10.0
        self._last_backpressure_log_ts: float = 0.0
        self.inflight_by_queue: Dict[str, int] = {}
        self.waiting_by_queue: Dict[str, List[Job]] = {}
        self.waiting_global: List[Job] = []
        # Optional on-disk spooling of waiting jobs to reduce memory pressure
        self._spool_enabled: bool = os.getenv("ROBIN_QUEUE_SPOOL", "1") != "0"
        self._spool_dir_source: str = "env"
        self.queue_spool_dir: Optional[str] = os.getenv("ROBIN_QUEUE_SPOOL_DIR")
        if not self.queue_spool_dir:
            self._spool_dir_source = "temp"
            self.queue_spool_dir = os.path.join(
                tempfile.gettempdir(), f"robin_queue_{os.getpid()}"
            )
        self._spool_counter: int = 0
        # per-sample tracking for GUI
        self.samples_by_id: Dict[str, Dict[str, Any]] = {}
        # Sample cleanup tracking to prevent unbounded memory growth
        self._sample_cleanup_interval: int = 5000  # Clean up every N jobs (conservative)
        self._jobs_since_last_cleanup: int = 0
        # shutdown flag for graceful termination
        self._shutdown_requested = False
        
        # CNV sample cache management
        self._cnv_cache_cleanup_interval: int = 1000  # Clean up CNV cache every N jobs
        self._jobs_since_cnv_cleanup: int = 0
        self._completed_samples: Set[str] = set()  # Track samples that have completed CNV analysis
        # defer setup/registration to async setup()

        # Preset controls how processing actors are created and their concurrency
        # None -> legacy per-type actors; otherwise use grouped Pool actors
        self.preset: Optional[str] = preset
        self.using_pools: bool = False

        # Reference genome for SNP calling and other analyses
        self.reference: Optional[str] = reference

        # Reference genome status (logged at INFO level)
        if self.reference:
            pass  # Reference genome available
        else:
            pass  # No reference genome

        # Keep references to any CPU limiter actors
        self._cpu_limiters: List[Any] = []
        
        # CSV buffering for performance (write in batches instead of per-job)
        self._csv_buffer: List[Dict[str, Any]] = []
        self._csv_buffer_size: int = 100  # Flush every 100 jobs
        self._csv_buffer_max_size: int = 500  # Maximum buffer size before forced flush
        # Lock created lazily to avoid serialization issues in Ray actors
        self._csv_buffer_lock: Optional[asyncio.Lock] = None
        
        # Stats caching for performance (avoid recomputing stats on every GUI update)
        self._stats_cache: Optional[Dict[str, Any]] = None
        self._stats_cache_time: float = 0.0
        self._stats_cache_ttl: float = 5.0  # Cache stats for 5 seconds
        
        # Batching support
        self.enable_batching = enable_batching
        self.sample_batcher = SampleJobBatcher() if enable_batching else None

        # Work directory (set by set_work_dir(); used for failed_jobs log)
        self.work_dir: Optional[str] = None
        # Bounded log of failed job details for tracking and OOM diagnosis (last 5000)
        self._failed_jobs_log: deque = deque(maxlen=5000)
        # Per-job timeout: cancel jobs running longer than this (0 = disabled)
        self.job_timeout_seconds: int = max(0, int(job_timeout_seconds))
    
    def _get_csv_buffer_lock(self) -> Optional[asyncio.Lock]:
        """Get or create the CSV buffer lock lazily."""
        if self._csv_buffer_lock is None and asyncio:
            self._csv_buffer_lock = asyncio.Lock()
        return self._csv_buffer_lock

    def get_reference(self) -> Optional[str]:
        """Get the reference genome path."""
        return self.reference

    def get_target_panel(self) -> str:
        """Get the target panel for fusion analysis."""
        return self.target_panel
    
    def set_work_dir(self, work_dir: Optional[str]) -> None:
        """Set the work directory for the coordinator."""
        self.work_dir = work_dir
        # If no explicit spool dir was provided, prefer work_dir for spooling
        if (
            self._spool_enabled
            and self._spool_dir_source == "temp"
            and work_dir
            and not os.getenv("ROBIN_QUEUE_SPOOL_DIR")
        ):
            self.queue_spool_dir = os.path.join(work_dir, ".robin_queue")

    def _write_failed_jobs_file(self) -> Optional[str]:
        """Write failed jobs log to work_dir/failed_jobs.csv for inspection. Returns path or None."""
        if not self._failed_jobs_log or not getattr(self, "work_dir", None) or not self.work_dir:
            return None
        try:
            import csv
            path = os.path.join(self.work_dir, "failed_jobs.csv")
            os.makedirs(self.work_dir, exist_ok=True)
            with open(path, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(["sample_id", "filepath", "job_type", "failure_kind", "is_oom", "error"])
                for e in self._failed_jobs_log:
                    w.writerow([
                        e.get("sample_id", ""),
                        e.get("filepath", ""),
                        e.get("job_type", ""),
                        e.get("failure_kind", ""),
                        "yes" if e.get("is_oom") else "no",
                        (e.get("error") or "").replace("\r", " ").replace("\n", " "),
                    ])
            return path
        except Exception as e:
            logging.warning(f"Could not write failed_jobs.csv: {e}")
            return None

    def _ensure_spool_dir(self) -> None:
        if not self._spool_enabled or not self.queue_spool_dir:
            return
        try:
            os.makedirs(self.queue_spool_dir, exist_ok=True)
        except Exception:
            pass

    def _spool_job(self, job: Job, sample_id: Optional[str]) -> Any:
        if not self._spool_enabled:
            return job
        self._ensure_spool_dir()
        if not self.queue_spool_dir:
            return job
        self._spool_counter += 1
        fname = f"job_{self._spool_counter}_{job.job_id}_{uuid.uuid4().hex}.pkl"
        path = os.path.join(self.queue_spool_dir, fname)
        try:
            with open(path, "wb") as fh:
                pickle.dump(job, fh, protocol=pickle.HIGHEST_PROTOCOL)
            return {"spool_path": path, "sample_id": sample_id or "unknown"}
        except Exception:
            return job

    def _unspool_job(self, entry: Any) -> Tuple[Optional[Job], Optional[str]]:
        if isinstance(entry, dict) and entry.get("spool_path"):
            path = entry.get("spool_path")
            sample_id = entry.get("sample_id") or "unknown"
            try:
                with open(path, "rb") as fh:
                    job = pickle.load(fh)
                try:
                    os.remove(path)
                except Exception:
                    pass
                return job, sample_id
            except Exception:
                return None, sample_id
        if isinstance(entry, Job):
            try:
                sid = entry.context.get_sample_id()
            except Exception:
                sid = "unknown"
            return entry, sid
        return None, None

    # ----- handler registration API -----
    async def register_handler(self, job_type: str, remote_func) -> None:
        opts = RESOURCE_HINTS.get(job_type, {})
        proc = self.processors.get(job_type)
        if proc is not None:
            await proc.update_handler.remote(remote_func, opts)

    async def setup(self) -> None:
        """Async setup to register default handlers and start drain loop."""
        # remote function handles (job_type -> remote function)
        registrations = [
            ("preprocessing", preprocessing_handler_remote),
            ("bed_conversion", bed_conversion_handler_remote),
            ("mgmt", mgmt_handler_remote),
            ("cnv", cnv_handler_remote),
            ("target", target_handler_remote),
            ("fusion", fusion_handler_remote),
            ("sturgeon", sturgeon_handler_remote),
            ("nanodx", nanodx_handler_remote),
            ("pannanodx", pannanodx_handler_remote),
            ("random_forest", random_forest_handler_remote),
            ("igv_bam", ray.remote(_wrap_real_handler(_igv_bam_handler, "igv_bam"))),
            (
                "snp_analysis",
                ray.remote(_wrap_real_handler(_snp_analysis_handler, "snp_analysis")),
            ),
            (
                "target_bam_finalize",
                ray.remote(
                    _wrap_real_handler(
                        _target_bam_finalize_handler, "target_bam_finalize"
                    )
                ),
            ),
        ]

        preset = (self.preset or "").lower().strip()
        if preset in {"p2i", "standard"}:
            # Optionally reserve CPUs to cap global availability
            try:
                total_cpus = float((ray.cluster_resources() or {}).get("CPU", 0))
            except Exception:
                total_cpus = 0.0
            desired_cap = 2.0 if preset == "p2i" else 6.0  # Increased from 4.0 to 6.0 for additional pools
            reserve = max(0.0, total_cpus - desired_cap)
            try:
                if reserve >= 1.0:
                    # Reserve floor(reserve) CPUs via a dummy actor
                    reserve_int = int(reserve)

                    @ray.remote
                    class _CpuLimiter:
                        async def ping(self):
                            return True

                    try:
                        limiter = _CpuLimiter.options(num_cpus=reserve_int).remote()
                    except Exception:
                        limiter = _CpuLimiter.remote()
                    # Keep reference so it's not GC'd
                    self._cpu_limiters.append(limiter)
            except Exception:
                pass

            # Grouped Pool actors - optimized resource allocation
            groups = {
                "prep": ["preprocessing"],  # Separate preprocessing for independence
                "bed_conversion": ["bed_conversion"],  # Separate bed_conversion for independence
                "cnv": ["cnv"],  # CNV gets its own CPU (most CPU-intensive)
                "analysis": ["mgmt", "target", "fusion"],  # Lightweight analysis types share a CPU
                "classif": ["sturgeon", "nanodx", "pannanodx"],  # Fast classification jobs
                "rf": ["random_forest"],  # Random forest needs its own actor (slow/blocking)
                "slow": [
                    "igv_bam",
                    "snp_analysis",
                    "target_bam_finalize",
                ],  # Slow pool for igv_bam and snp_analysis jobs
            }

            # Determine concurrency per pool based on preset
            def pool_parallel(name: str) -> int:
                if preset == "p2i":
                    # Concurrency 1 everywhere
                    return 1
                # standard preset - optimized resource allocation
                if name == "cnv":
                    # CNV gets dedicated concurrency (most CPU-intensive)
                    return max(1, int(self.analysis_workers))
                if name == "analysis":
                    # Lightweight analysis types share concurrency
                    return max(1, int(self.analysis_workers))
                if name == "classif":
                    # Fast classification jobs get moderate concurrency
                    return max(1, int(self.analysis_workers))
                if name == "rf":
                    # Random forest gets its own concurrency (slow/blocking)
                    return max(1, int(self.analysis_workers))
                # Give preprocessing and bed_conversion pools moderate concurrency
                if name in {"prep", "bed_conversion"}:
                    return max(1, int(self.analysis_workers))
                return 1

            # Create pools and register handlers
            pools: Dict[str, Any] = {}
            for pool_name, job_types in groups.items():
                par = pool_parallel(pool_name)
                actor_name = f"pool_{pool_name}"
                job_timeout = getattr(self, "job_timeout_seconds", 0) or 0
                try:
                    pool = Pool.options(
                        name=actor_name, max_concurrency=max(1, int(par)), num_cpus=0
                    ).remote(pool_name, par, job_timeout)
                except Exception:
                    try:
                        # Some Ray versions don't allow num_cpus=0; use a tiny fractional CPU to avoid reservation deadlock
                        pool = Pool.options(
                            name=actor_name,
                            max_concurrency=max(1, int(par)),
                            num_cpus=0.001,
                        ).remote(pool_name, par, job_timeout)
                    except Exception:
                        pool = Pool.options(max_concurrency=max(1, int(par))).remote(
                            pool_name, par, job_timeout
                        )
                pools[pool_name] = pool
                # Wire coordinator callback
                try:
                    # Fire-and-forget to avoid blocking if actor cannot start immediately
                    pool.set_coordinator_name.remote("robin_coordinator")
                except Exception:
                    pass
                # Register handlers on the pool and map processors
                for jt, rf in registrations:
                    if jt in job_types:
                        opts = RESOURCE_HINTS.get(jt, {})
                        try:
                            await pool.register_handler.remote(jt, rf, opts)
                        except Exception:
                            pass
                        self.processors[jt] = pool
            self.using_pools = True
        elif preset == "high":
            # High-performance: per-type actors (legacy model), actors reserve 1 CPU by default
            for jt, rf in registrations:
                opts = RESOURCE_HINTS.get(jt, {})
                max_conc = 1
                if jt in {"mgmt", "cnv", "target", "fusion"}:
                    max_conc = max(1, int(self.analysis_workers))
                name = f"typeproc_{jt}"
                try:
                    proc = TypeProcessor.options(
                        name=name, max_concurrency=max_conc
                    ).remote(jt, rf, opts)
                except Exception:
                    proc = TypeProcessor.options(max_concurrency=max_conc).remote(
                        jt, rf, opts
                    )
                self.processors[jt] = proc
            self.using_pools = False
        else:
            # Legacy: one TypeProcessor per job type
            for jt, rf in registrations:
                opts = RESOURCE_HINTS.get(jt, {})
                # Set actor concurrency via .options on the actor itself
                max_conc = 1
                if jt in {"mgmt", "cnv", "target", "fusion"}:
                    max_conc = max(1, int(self.analysis_workers))
                name = f"typeproc_{jt}"
                try:
                    proc = TypeProcessor.options(
                        name=name, max_concurrency=max_conc
                    ).remote(jt, rf, opts)
                except Exception:
                    proc = TypeProcessor.options(max_concurrency=max_conc).remote(
                        jt, rf, opts
                    )
                self.processors[jt] = proc
            self.using_pools = False
        # start drain loop to observe completions
        try:
            asyncio.create_task(self._drain_loop())
        except Exception:
            pass
        
        # Start batch timeout loop if batching is enabled
        if self.enable_batching and self.sample_batcher:
            try:
                asyncio.create_task(self._batch_timeout_loop())
            except Exception:
                pass

    def _total_inflight(self) -> int:
        try:
            return len(self.active)
        except Exception:
            return 0

    def _total_waiting(self) -> int:
        try:
            waiting_types = sum(len(v) for v in self.waiting_by_type_global.values())
        except Exception:
            waiting_types = 0
        try:
            waiting_queues = sum(len(v) for v in self.waiting_by_queue.values())
        except Exception:
            waiting_queues = 0
        try:
            waiting_global = len(self.waiting_global)
        except Exception:
            waiting_global = 0
        return waiting_types + waiting_queues + waiting_global

    def _queue_inflight_cap(self, queue_name: str) -> int:
        return int(self.max_inflight_per_queue.get(queue_name, self.max_inflight_per_type))

    def _queue_waiting_cap(self, queue_name: str) -> int:
        return int(self.max_waiting_per_queue.get(queue_name, self.max_inflight_per_type * 20))

    async def _wait_for_global_capacity(self) -> None:
        while self._total_waiting() >= self.max_total_waiting:
            await asyncio.sleep(0.05)

    async def _wait_for_queue_capacity(self, queue_name: str) -> None:
        while len(self.waiting_by_queue.get(queue_name, [])) >= self._queue_waiting_cap(queue_name):
            await asyncio.sleep(0.05)

    async def can_accept_jobs(self, n: int = 1) -> bool:
        try:
            n = int(n or 0)
        except Exception:
            n = 0
        return (self._total_waiting() + n) < self.max_total_waiting

    async def _dispatch_ready_job(self, job: Job, sample_id: str, from_waiting: bool = False) -> None:
        # Ensure global waiting does not grow without bound for new submissions
        if not from_waiting:
            await self._wait_for_global_capacity()

        # Global inflight cap
        if self._total_inflight() >= self.max_total_inflight:
            if self._debug_backpressure:
                now = time.time()
                if (now - self._last_backpressure_log_ts) >= self._debug_backpressure_interval_s:
                    self._last_backpressure_log_ts = now
                    logging.info(
                        "[ROBIN][backpressure] global inflight cap hit: "
                        f"inflight={self._total_inflight()}/{self.max_total_inflight}, "
                        f"waiting_global={len(self.waiting_global)}, "
                        f"waiting_serialized={self._total_waiting()}"
                    )
            self.waiting_global.append(self._spool_job(job, sample_id))
            try:
                sid_pending = sample_id or "unknown"
                if sid_pending != "unknown":
                    ent_pending = self.samples_by_id.get(sid_pending)
                    if ent_pending is None:
                        ent_pending = {
                            "sample_id": sid_pending,
                            "active_jobs": 0,
                            "pending_jobs": 0,
                            "total_jobs": 0,
                            "completed_jobs": 0,
                            "failed_jobs": 0,
                            "job_types": set(),
                            "last_seen": time.time(),
                        }
                        self.samples_by_id[sid_pending] = ent_pending
                    ent_pending["pending_jobs"] = ent_pending.get("pending_jobs", 0) + 1
                    ent_pending["last_seen"] = time.time()
            except Exception:
                pass
            return

        q = job_queue_of(job.job_type)
        # Per-queue inflight cap
        if int(self.inflight_by_queue.get(q, 0)) >= self._queue_inflight_cap(q):
            await self._wait_for_queue_capacity(q)
            if self._debug_backpressure:
                now = time.time()
                if (now - self._last_backpressure_log_ts) >= self._debug_backpressure_interval_s:
                    self._last_backpressure_log_ts = now
                    logging.info(
                        "[ROBIN][backpressure] queue inflight cap hit: "
                        f"queue={q}, inflight={self.inflight_by_queue.get(q,0)}/{self._queue_inflight_cap(q)}, "
                        f"waiting_queue={len(self.waiting_by_queue.get(q, []))}/{self._queue_waiting_cap(q)}"
                    )
            self.waiting_by_queue.setdefault(q, []).append(
                self._spool_job(job, sample_id)
            )
            try:
                sid_pending = sample_id or "unknown"
                if sid_pending != "unknown":
                    ent_pending = self.samples_by_id.get(sid_pending)
                    if ent_pending is None:
                        ent_pending = {
                            "sample_id": sid_pending,
                            "active_jobs": 0,
                            "pending_jobs": 0,
                            "total_jobs": 0,
                            "completed_jobs": 0,
                            "failed_jobs": 0,
                            "job_types": set(),
                            "last_seen": time.time(),
                        }
                        self.samples_by_id[sid_pending] = ent_pending
                    ent_pending["pending_jobs"] = ent_pending.get("pending_jobs", 0) + 1
                    ent_pending["last_seen"] = time.time()
            except Exception:
                pass
            return

        # Per-type inflight cap
        inflight_for_type = int(self.inflight_by_type.get(job.job_type, 0))
        if inflight_for_type >= self.max_inflight_per_type:
            if self._debug_backpressure:
                now = time.time()
                if (now - self._last_backpressure_log_ts) >= self._debug_backpressure_interval_s:
                    self._last_backpressure_log_ts = now
                    logging.info(
                        "[ROBIN][backpressure] per-type inflight cap hit: "
                        f"type={job.job_type}, inflight={inflight_for_type}/{self.max_inflight_per_type}, "
                        f"waiting_global={len(self.waiting_by_type_global.get(job.job_type, []))}"
                    )
            self.waiting_by_type_global.setdefault(job.job_type, []).append(
                self._spool_job(job, sample_id)
            )
            try:
                sid_pending = sample_id or "unknown"
                if sid_pending != "unknown":
                    ent_pending = self.samples_by_id.get(sid_pending)
                    if ent_pending is None:
                        ent_pending = {
                            "sample_id": sid_pending,
                            "active_jobs": 0,
                            "pending_jobs": 0,
                            "total_jobs": 0,
                            "completed_jobs": 0,
                            "failed_jobs": 0,
                            "job_types": set(),
                            "last_seen": time.time(),
                        }
                        self.samples_by_id[sid_pending] = ent_pending
                    ent_pending["pending_jobs"] = ent_pending.get("pending_jobs", 0) + 1
                    ent_pending["last_seen"] = time.time()
            except Exception:
                pass
            return

        # Found processor for job type
        proc = self.processors.get(job.job_type)
        if proc is None:
            return  # Skip jobs without processors

        # Mark submission only when actually dispatching to a processor
        self.total_enqueued += 1
        try:
            self.submitted_by_type[job.job_type] = (
                self.submitted_by_type.get(job.job_type, 0) + 1
            )
        except Exception:
            pass
        self.active[job.job_id] = {
            "job_type": job.job_type,
            "filepath": job.context.filepath,
            "queue": q,
            "start_time": time.time(),
        }
        # Update per-sample totals/active only when a real sample_id is known
        try:
            sid = sample_id or "unknown"
            if sid != "unknown":
                ent = self.samples_by_id.get(sid)
                if ent is None:
                    ent = {
                        "sample_id": sid,
                        "active_jobs": 0,
                        "pending_jobs": 0,
                        "total_jobs": 0,
                        "completed_jobs": 0,
                        "failed_jobs": 0,
                        "job_types": set(),
                        "last_seen": time.time(),
                    }
                    self.samples_by_id[sid] = ent
                ent["total_jobs"] += 1
                ent["active_jobs"] += 1
                try:
                    ent["job_types"].add(job.job_type)
                except Exception:
                    pass
                ent["last_seen"] = time.time()
        except Exception:
            pass
        if getattr(self, "using_pools", False):
            try:
                proc.enqueue.remote(job)
            except Exception:
                pass
        else:
            ref = proc.process.remote(job)
            self._inflight[ref] = job
        self.inflight_by_type[job.job_type] = inflight_for_type + 1
        self.inflight_by_queue[q] = int(self.inflight_by_queue.get(q, 0)) + 1

    async def _drain_waiting_jobs(self, queue_name: Optional[str], job_type: str) -> None:
        # Try one global waiting job first to avoid starvation
        if self.waiting_global and self._total_inflight() < self.max_total_inflight:
            entry = self.waiting_global.pop(0)
            nxt_global, sid = self._unspool_job(entry)
            if nxt_global is None:
                return
            try:
                if sid != "unknown":
                    ent_pending = self.samples_by_id.get(sid)
                    if ent_pending and ent_pending.get("pending_jobs", 0) > 0:
                        ent_pending["pending_jobs"] -= 1
                        ent_pending["last_seen"] = time.time()
            except Exception:
                pass
            await self._dispatch_ready_job(nxt_global, sid, from_waiting=True)

        # Then try one queued-by-queue job
        if queue_name:
            queue_jobs = self.waiting_by_queue.get(queue_name, [])
            if queue_jobs and int(self.inflight_by_queue.get(queue_name, 0)) < self._queue_inflight_cap(queue_name):
                entry = queue_jobs.pop(0)
                if not queue_jobs:
                    self.waiting_by_queue.pop(queue_name, None)
                nxt_queue, sid = self._unspool_job(entry)
                if nxt_queue is None:
                    return
                try:
                    if sid != "unknown":
                        ent_pending = self.samples_by_id.get(sid)
                        if ent_pending and ent_pending.get("pending_jobs", 0) > 0:
                            ent_pending["pending_jobs"] -= 1
                            ent_pending["last_seen"] = time.time()
                except Exception:
                    pass
                await self._dispatch_ready_job(nxt_queue, sid, from_waiting=True)

        # Finally, release one per-type waiting job
        queue_global = self.waiting_by_type_global.get(job_type, [])
        if queue_global and int(self.inflight_by_type.get(job_type, 0)) < self.max_inflight_per_type:
            entry = queue_global.pop(0)
            if not queue_global:
                self.waiting_by_type_global.pop(job_type, None)
            nxt_job_g, sid = self._unspool_job(entry)
            if nxt_job_g is None:
                return
            try:
                if sid != "unknown":
                    ent_pending = self.samples_by_id.get(sid)
                    if ent_pending and ent_pending.get("pending_jobs", 0) > 0:
                        ent_pending["pending_jobs"] -= 1
                        ent_pending["last_seen"] = time.time()
            except Exception:
                pass
            await self._dispatch_ready_job(nxt_job_g, sid, from_waiting=True)

    # ----- submission & lifecycle -----
    async def submit_jobs(self, jobs: List[Job]) -> None:
        # Check if shutdown has been requested
        if self._shutdown_requested:
            return
        
        # Invalidate stats cache since we're submitting new jobs
        self._stats_cache = None
        
        # Process jobs through batcher if batching is enabled
        if self.enable_batching and self.sample_batcher:
            # Separate preprocessing jobs from other jobs
            preprocessing_jobs = [job for job in jobs if job.job_type == "preprocessing"]
            other_jobs = [job for job in jobs if job.job_type != "preprocessing"]
            
            # Submit preprocessing jobs directly (no batching)
            if preprocessing_jobs:
                await self._submit_jobs_internal(preprocessing_jobs)
            
            # Process other jobs through batcher
            for job in other_jobs:
                try:
                    sid_pending = job.context.get_sample_id() if job.context else "unknown"
                except Exception:
                    sid_pending = "unknown"
                try:
                    if sid_pending != "unknown":
                        ent_pending = self.samples_by_id.get(sid_pending)
                        if ent_pending is None:
                            ent_pending = {
                                "sample_id": sid_pending,
                                "active_jobs": 0,
                                "pending_jobs": 0,
                                "total_jobs": 0,
                                "completed_jobs": 0,
                                "failed_jobs": 0,
                                "job_types": set(),
                                "last_seen": time.time(),
                            }
                            self.samples_by_id[sid_pending] = ent_pending
                        ent_pending["pending_jobs"] = ent_pending.get("pending_jobs", 0) + 1
                        ent_pending["last_seen"] = time.time()
                except Exception:
                    pass
                batches = await self.sample_batcher.add_job(job)
                if batches:
                    try:
                        for batch in batches:
                            sid_batch = batch.sample_id
                            if sid_batch != "unknown":
                                ent_pending = self.samples_by_id.get(sid_batch)
                                if ent_pending and ent_pending.get("pending_jobs", 0) > 0:
                                    ent_pending["pending_jobs"] = max(
                                        0,
                                        ent_pending.get("pending_jobs", 0) - len(batch.contexts),
                                    )
                                    ent_pending["last_seen"] = time.time()
                    except Exception:
                        pass
                    # Convert BatchedJob to regular Job for processing
                    regular_jobs = self._convert_batched_jobs(batches)
                    await self._submit_jobs_internal(regular_jobs)
                # If no batches were created, the job is waiting in the batcher
        else:
            await self._submit_jobs_internal(jobs)
    
    async def _submit_jobs_internal(self, jobs: List[Job]) -> None:
        """Internal job submission logic"""
        for job in jobs:
            sample_id = job.context.get_sample_id()
            if job.job_type in DEDUP_TYPES:
                key = (sample_id, job.job_type)
                run = self.running.get(key, 0)
                pend = self.pending.get(key, 0)
                if run + pend >= 2 or pend >= 1:
                    self.total_skipped += 1
                    continue
                self.pending[key] = pend + 1
            # Per-sample per-type serialization; queue instead of overlapping
            if job.job_type in SERIALIZE_BY_TYPE_PER_SAMPLE:
                sid = sample_id or "unknown"
                key_ts = (job.job_type, sid)
                ra = self.running_by_type_sample.get(key_ts, 0)
                pa = self.pending_by_type_sample.get(key_ts, 0)
                if (ra + pa) > 0:
                    self.waiting_by_type_sample.setdefault(key_ts, []).append(job)
                    self.pending_by_type_sample[key_ts] = pa + 1
                    try:
                        if sid != "unknown":
                            ent_pending = self.samples_by_id.get(sid)
                            if ent_pending is None:
                                ent_pending = {
                                    "sample_id": sid,
                                    "active_jobs": 0,
                                    "pending_jobs": 0,
                                    "total_jobs": 0,
                                    "completed_jobs": 0,
                                    "failed_jobs": 0,
                                    "job_types": set(),
                                    "last_seen": time.time(),
                                }
                                self.samples_by_id[sid] = ent_pending
                            ent_pending["pending_jobs"] = ent_pending.get("pending_jobs", 0) + 1
                            ent_pending["last_seen"] = time.time()
                    except Exception:
                        pass
                    continue
                else:
                    self.running_by_type_sample[key_ts] = ra + 1

            await self._dispatch_ready_job(job, sample_id)

    async def submit_sample_job(
        self,
        sample_dir: str,
        job_type: str,
        sample_id: str = None,
        force_regenerate: bool = False,
    ) -> bool:
        """
        Submit a job for an existing sample directory.

        This allows users to manually trigger specific job types for samples
        that have already been processed or need reprocessing.

        Args:
            sample_dir: Path to the sample directory
            job_type: Type of job to run (e.g., 'igv_bam')
            sample_id: Optional sample ID (defaults to directory name)
            force_regenerate: If True, force regeneration even if output exists

        Returns:
            True if job was successfully submitted, False otherwise
        """
        try:
            if sample_id is None:
                sample_id = Path(sample_dir).name

            # Create a context for this sample
            context = WorkflowContext(
                filepath=sample_dir,
                metadata={
                    "sample_id": sample_id,
                    "sample_dir": sample_dir,
                    "work_dir": os.path.dirname(
                        sample_dir
                    ),  # Set the work_dir to parent directory
                    "force_regenerate": force_regenerate,  # Add force_regenerate flag
                    "bam_metadata": {"sample_id": sample_id},
                },
            )

            # Create a job
            job = Job(
                job_id=next(_job_id_counter),
                job_type=job_type,
                origin="manual",  # Use 'origin' instead of 'queue'
                workflow=[f"slow:{job_type}"],
                step=0,
                context=context,
            )

            # Submit the job
            await self.submit_jobs([job])

            return True

        except Exception:
            return False

    async def submit_snp_analysis_job(
        self,
        sample_dir: str,
        sample_id: str = None,
        reference: str = None,
        threads: int = 4,
        force_regenerate: bool = False,
    ) -> bool:
        """
        Submit a SNP analysis job for an existing sample directory.

        This is a convenience method specifically for SNP analysis that ensures
        all required metadata is properly set up.

        Args:
            sample_dir: Path to the sample directory
            sample_id: Optional sample ID (defaults to directory name)
            reference: Path to reference genome (optional, will auto-detect if not provided)
            threads: Number of threads to use for processing (default: 4)
            force_regenerate: Whether to force regeneration of existing results (default: False)

        Returns:
            True if job was successfully submitted, False otherwise
        """
        logger = logging.getLogger("robin.workflow")
        try:
            if sample_id is None:
                sample_id = Path(sample_dir).name
            logger.info(
                "Submitting SNP analysis job: sample_id=%s sample_dir=%s reference=%s threads=%s force_regenerate=%s",
                sample_id,
                sample_dir,
                reference,
                threads,
                force_regenerate,
            )

            # Determine target panel from master.csv if present
            target_panel = None
            try:
                master_csv = Path(sample_dir) / "master.csv"
                if master_csv.exists():
                    import csv

                    with master_csv.open("r", newline="") as fh:
                        reader = csv.DictReader(fh)
                        first_row = next(reader, None)
                        if first_row:
                            panel = first_row.get("analysis_panel", "").strip()
                            if panel:
                                target_panel = panel
            except Exception:
                pass

            # Create a context for this sample with SNP-specific metadata
            metadata = {
                "sample_id": sample_id,
                "sample_dir": sample_dir,
                "work_dir": os.path.dirname(
                    sample_dir
                ),  # Set the work_dir to parent directory
                "threads": threads,
                "force_regenerate": force_regenerate,
                "bam_metadata": {"sample_id": sample_id},
            }

            if target_panel:
                metadata["target_panel"] = target_panel

            # Add reference genome if provided
            if reference:
                metadata["reference"] = reference
            else:
                logger.warning(
                    "SNP submission for sample_id=%s has no reference in metadata; downstream handler may skip.",
                    sample_id,
                )

            context = WorkflowContext(filepath=sample_dir, metadata=metadata)

            # Create a job
            job = Job(
                job_id=next(_job_id_counter),
                job_type="snp_analysis",
                origin="manual",  # Use 'origin' instead of 'queue'
                workflow=["slow:snp_analysis"],
                step=0,
                context=context,
            )

            # Submit the job
            await self.submit_jobs([job])
            logger.info("SNP analysis job queued: sample_id=%s job_id=%s", sample_id, job.job_id)

            return True

        except Exception as e:
            logger.error(
                "Failed to submit SNP analysis job: sample_id=%s sample_dir=%s error=%s",
                sample_id,
                sample_dir,
                e,
                exc_info=True,
            )
            return False

    async def submit_target_bam_finalize_job(
        self,
        sample_dir: str,
        sample_id: str = None,
        target_panel: str = None,
        reference: str = None,
    ) -> bool:
        """
        Submit a target BAM finalization job for an existing sample directory.
        """
        logger = logging.getLogger("robin.workflow")
        try:
            if sample_id is None:
                sample_id = Path(sample_dir).name
            logger.info(
                "Submitting target_bam_finalize job: sample_id=%s sample_dir=%s target_panel=%s",
                sample_id,
                sample_dir,
                target_panel,
            )

            if target_panel is None:
                try:
                    import csv
                    master_csv = Path(sample_dir) / "master.csv"
                    if master_csv.exists():
                        with master_csv.open("r", newline="") as fh:
                            reader = csv.DictReader(fh)
                            first_row = next(reader, None)
                            if first_row:
                                panel = first_row.get("analysis_panel", "").strip()
                                if panel:
                                    target_panel = panel
                except Exception:
                    pass

            metadata = {
                "sample_id": sample_id,
                "sample_dir": sample_dir,
                "work_dir": os.path.dirname(sample_dir),
                "bam_metadata": {"sample_id": sample_id},
            }

            if target_panel:
                metadata["target_panel"] = target_panel
            if reference:
                metadata["reference"] = reference

            context = WorkflowContext(filepath=sample_dir, metadata=metadata)

            job = Job(
                job_id=next(_job_id_counter),
                job_type="target_bam_finalize",
                origin="manual",
                workflow=["slow:target_bam_finalize"],
                step=0,
                context=context,
            )

            await self.submit_jobs([job])
            logger.info(
                "target_bam_finalize job queued: sample_id=%s job_id=%s",
                sample_id,
                job.job_id,
            )
            return True

        except Exception as e:
            logger.error(
                "Failed to submit target_bam_finalize job: sample_id=%s sample_dir=%s error=%s",
                sample_id,
                sample_dir,
                e,
                exc_info=True,
            )
            return False

    def is_sample_ready_for_snp_analysis(
        self, sample_dir: str
    ) -> tuple[bool, list[str]]:
        """
        Check if a sample directory is ready for SNP analysis.

        Args:
            sample_dir: Path to the sample directory

        Returns:
            Tuple of (is_ready, missing_files) where is_ready is a boolean
            and missing_files is a list of missing required files
        """
        required_files = ["target.bam", "targets_exceeding_threshold.bed"]

        missing_files = []
        for filename in required_files:
            file_path = os.path.join(sample_dir, filename)
            if not os.path.exists(file_path):
                missing_files.append(filename)

        is_ready = len(missing_files) == 0
        return is_ready, missing_files

    async def _cleanup_cnv_cache(self) -> None:
        """
        Clean up CNV sample cache for completed samples to prevent unbounded memory growth.
        """
        if clear_sample_cache is None:
            return  # CNV module not available
        
        try:
            # Get current cache stats
            cache_stats = get_sample_cache_stats() if get_sample_cache_stats else {}
            current_time = time.time()
            
            # Clean up cache for samples that haven't been accessed recently (10 minutes)
            samples_to_cleanup = []
            for sample_id, details in cache_stats.get('sample_details', {}).items():
                last_accessed = details.get('last_accessed', 0)
                if current_time - last_accessed > 600:  # 10 minutes
                    samples_to_cleanup.append(sample_id)
            
            # Limit cleanup to 5 samples at a time to avoid disruption
            samples_to_cleanup = samples_to_cleanup[:5]
            
            # Clean up old samples
            for sample_id in samples_to_cleanup:
                try:
                    clear_sample_cache(sample_id)
                except Exception:
                    pass
            
            if samples_to_cleanup:
                pass  # Cache cleanup completed
                
        except Exception:
            pass  # Silently continue on cache cleanup errors

    async def _cleanup_completed_samples(self) -> None:
        """
        Clean up samples that have completed all their jobs.
        This prevents unbounded memory growth in samples_by_id dictionary.
        
        Conservative cleanup: only removes samples that have been idle for 5+ minutes.
        """
        if not self.samples_by_id:
            return
        
        # Only cleanup if we have a large number of samples (> 1000)
        if len(self.samples_by_id) < 1000:
            return
        
        # Identify samples with no active jobs and that haven't been seen recently
        current_time = time.time()
        samples_to_remove = []
        
        for sid, ent in self.samples_by_id.items():
            # Keep samples with active jobs
            if ent.get("active_jobs", 0) > 0:
                continue
            
            # Remove samples that completed all jobs and haven't been active for 5 minutes
            # This is very conservative to avoid removing data that might still be needed
            last_seen = ent.get("last_seen", 0)
            if current_time - last_seen > 300:  # 5 minutes
                samples_to_remove.append(sid)
        
        # Limit cleanup to 100 samples at a time to avoid disruption
        samples_to_remove = samples_to_remove[:100]
        
        # Remove completed samples
        for sid in samples_to_remove:
            self.samples_by_id.pop(sid, None)
    
    async def _finalize_target_bams(self) -> None:
        """
        Finalize target BAM accumulation for all samples that have batch BAMs pending merge.
        This should be called during shutdown to ensure all target.bam files are created.
        Uses Ray tasks to offload the blocking I/O work to separate workers.
        """
        # Try to get work_dir from Coordinator attribute first
        work_dir = getattr(self, 'work_dir', None)
        
        # If not set, try to extract from active jobs' metadata
        if not work_dir:
            for job in self._inflight.values():
                try:
                    work_dir = job.context.metadata.get("work_dir")
                    if work_dir:
                        break
                except Exception:
                    continue
        
        # If still not found, try to extract from samples_by_id if they have work_dir stored
        if not work_dir:
            for sample_data in self.samples_by_id.values():
                try:
                    work_dir = sample_data.get("work_dir")
                    if work_dir:
                        break
                except Exception:
                    continue
        
        if not work_dir:
            print("[SHUTDOWN] No work directory found - skipping target BAM finalization")
            print("[SHUTDOWN] Work directory may not be set on Coordinator or available in job metadata")
            return
        
        try:
            # Create a Ray remote function for the blocking work
            @ray.remote
            def finalize_sample_task(sample_id: str, work_dir: str, target_panel: str):
                """Ray task wrapper for finalize_accumulation_for_sample."""
                from robin.analysis.target_analysis import finalize_accumulation_for_sample
                return finalize_accumulation_for_sample(
                    sample_id=sample_id,
                    work_dir=work_dir,
                    target_panel=target_panel
                )
            
            # Get all samples that have been processed
            samples_to_finalize = set()
            for sid in self.samples_by_id.keys():
                if sid and sid != "unknown":
                    samples_to_finalize.add(sid)
            
            if not samples_to_finalize:
                print("[SHUTDOWN] No samples to finalize")
                return
            
            total = len(samples_to_finalize)
            start_time = time.time()
            heartbeat_seconds = 10
            last_heartbeat = 0.0

            print(
                f"[SHUTDOWN] Finalizing target.bam for {total} sample(s) using work_dir: {work_dir}",
                flush=True,
            )
            print(
                "[SHUTDOWN] Progress updates will appear as each sample finalization completes.",
                flush=True,
            )
            
            # Submit all finalization tasks to Ray workers in parallel
            tasks = []
            max_queue_print = 10
            for i, sample_id in enumerate(samples_to_finalize):
                if i < max_queue_print:
                    print(f"[SHUTDOWN] Queued finalization for sample: {sample_id}", flush=True)
                elif i == max_queue_print:
                    remaining_to_queue = total - max_queue_print
                    if remaining_to_queue > 0:
                        print(f"[SHUTDOWN] Queued additional {remaining_to_queue} sample(s)...", flush=True)
                task_ref = finalize_sample_task.remote(
                    sample_id=sample_id,
                    work_dir=work_dir,
                    target_panel=self.target_panel
                )
                tasks.append((sample_id, task_ref))

            # Wait for tasks to complete with periodic progress output.
            # (Avoids an apparently "stuck" period while asyncio.gather() waits silently.)
            completed = 0
            task_ref_to_sample_id = {task_ref: sample_id for sample_id, task_ref in tasks}
            remaining_refs = [task_ref for _, task_ref in tasks]

            while remaining_refs:
                # Wait briefly for at least one task to become ready.
                try:
                    ready_refs, remaining_refs = await asyncio.to_thread(
                        ray.wait,
                        remaining_refs,
                        1,
                        5.0,  # timeout seconds; keeps the shutdown console responsive
                    )
                except Exception as e:
                    print(f"[SHUTDOWN] Warning: Error while waiting for target.bam finalization tasks: {e}", flush=True)
                    logging.warning(f"Error while waiting for target.bam finalization tasks: {e}")
                    break

                if ready_refs:
                    for ref in ready_refs:
                        sample_id = task_ref_to_sample_id.get(ref, "unknown")
                        try:
                            result = await asyncio.to_thread(ray.get, ref)
                            completed += 1

                            if isinstance(result, Exception):
                                print(
                                    f"[SHUTDOWN] Warning: Failed to finalize target.bam for {sample_id}: {result}",
                                    flush=True,
                                )
                                logging.warning(f"Failed to finalize target.bam for {sample_id}: {result}")
                                continue

                            if isinstance(result, dict) and result.get("status") == "error":
                                print(
                                    f"[SHUTDOWN] Warning: Finalization error for {sample_id}: {result.get('error', 'Unknown')}",
                                    flush=True,
                                )
                                logging.warning(
                                    f"Finalization error for {sample_id}: {result.get('error', 'Unknown')}"
                                )
                            else:
                                batch_count = 0
                                if isinstance(result, dict):
                                    batch_count = result.get("batch_files_merged", 0) or 0
                                if batch_count > 0:
                                    print(
                                        f"[SHUTDOWN] Completed: {sample_id} (merged {batch_count} batch files)",
                                        flush=True,
                                    )
                                else:
                                    print(f"[SHUTDOWN] Completed: {sample_id}", flush=True)
                        except Exception as e:
                            print(
                                f"[SHUTDOWN] Warning: Failed to finalize target.bam for {sample_id}: {e}",
                                flush=True,
                            )
                            logging.warning(f"Failed to finalize target.bam for {sample_id}: {e}")

                # Heartbeat: user-friendly output during long waits.
                now = time.time()
                if now - last_heartbeat >= heartbeat_seconds:
                    elapsed_seconds = int(now - start_time)
                    remaining = len(remaining_refs)
                    print(
                        f"[SHUTDOWN] Still finalizing target.bam... {completed}/{total} completed, {remaining} remaining (elapsed {elapsed_seconds}s)",
                        flush=True,
                    )
                    last_heartbeat = now

            print(f"[SHUTDOWN] Target BAM finalization: {completed}/{total} samples completed", flush=True)
        except Exception as e:
            print(f"[SHUTDOWN] Error during target BAM finalization: {e}")
            logging.warning(f"Error during target BAM finalization: {e}")
    
    async def _flush_csv_buffer(self) -> None:
        """Flush buffered CSV entries to disk."""
        if not self._csv_buffer:
            return
        
        # Use lock to prevent concurrent flushes
        lock = self._get_csv_buffer_lock()
        if lock:
            async with lock:
                buffer_to_write = self._csv_buffer[:]
                self._csv_buffer.clear()
        else:
            buffer_to_write = self._csv_buffer[:]
            self._csv_buffer.clear()
        
        if not buffer_to_write:
            return
        
        # Group entries by work_dir to minimize file operations
        entries_by_dir: Dict[str, List[Dict[str, Any]]] = {}
        for entry in buffer_to_write:
            work_dir = entry.get('work_dir', '')
            if work_dir not in entries_by_dir:
                entries_by_dir[work_dir] = []
            entries_by_dir[work_dir].append(entry)
        
        # Write each group to its respective file
        for work_dir, entries in entries_by_dir.items():
            try:
                out_path = (
                    os.path.join(work_dir, "job_durations.csv")
                    if work_dir
                    else "job_durations.csv"
                )
                
                # Create parent directory if needed
                try:
                    parent = os.path.dirname(out_path)
                    if parent:
                        os.makedirs(parent, exist_ok=True)
                except Exception:
                    pass
                
                # Check if we need to write header
                new_file = not os.path.exists(out_path)
                
                with open(out_path, "a", encoding="utf-8") as fh:
                    if new_file:
                        fh.write(
                            "timestamp_iso,sample_id,job_id,job_type,queue,filepath,duration_seconds,status,error\n"
                        )
                    
                    for entry in entries:
                        fh.write(
                            f"{entry['timestamp_iso']},{entry['sample_id']},{entry['job_id']},"
                            f"{entry['job_type']},{entry['queue']},{entry['filepath']},"
                            f"{entry['duration_seconds']},{entry['status']},{entry['error']}\n"
                        )
            except Exception:
                pass  # Silently continue on write errors

    async def _on_finish(
        self,
        job: Job,
        ok: bool,
        ctx: WorkflowContext,
        err: Optional[str],
        failure_kind: Optional[str] = None,
    ):
        # Infer structured failure kind when not provided (for filtering/alerting)
        if not ok and failure_kind is None and err:
            failure_kind = (
                "timed_out"
                if "timed out" in (err or "").lower()
                else "exception"
            )
        elif ok:
            failure_kind = None
        # capture start info for duration logging, then update dedup maps
        active_info = self.active.pop(job.job_id, None)
        start_time = None
        queue_name = None
        filepath = None
        try:
            if isinstance(active_info, dict):
                start_time = active_info.get("start_time")
                queue_name = active_info.get("queue")
                filepath = active_info.get("filepath")
        except Exception:
            pass
        # update dedup maps
        if job.job_type in DEDUP_TYPES:
            key = (job.context.get_sample_id(), job.job_type)
            self.running[key] = max(0, self.running.get(key, 0) - 1)
            if self.running[key] == 0:
                self.running.pop(key, None)

        # update stats
        # Detect deliberate skip (large BAM)
        is_skipped = False
        try:
            is_skipped = (job.job_type == "preprocessing") and bool(
                ctx.metadata.get("skip_downstream", False)
            )
        except Exception:
            is_skipped = False
        status = (
            "completed"
            if (ok and not is_skipped)
            else ("skipped" if is_skipped else "failed")
        )

        # Append per-job duration entry to buffer (for batched writes)
        try:
            end_time = time.time()
            duration = (end_time - float(start_time)) if start_time else None
            sample_id = None
            try:
                sample_id = (
                    ctx.get_sample_id() if hasattr(ctx, "get_sample_id") else "unknown"
                )
            except Exception:
                sample_id = "unknown"
            work_dir = None
            try:
                work_dir = ctx.metadata.get("work_dir") or job.context.metadata.get(
                    "work_dir"
                )
            except Exception:
                work_dir = None
            
            ts_iso = time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(end_time))
            dur_str = (
                f"{duration:.3f}" if isinstance(duration, (int, float)) else ""
            )
            fp = filepath or (ctx.filepath if hasattr(ctx, "filepath") else "")
            err_short = (err or "").replace("\n", " ").replace(",", " ")[:500]
            
            # Add entry to buffer
            self._csv_buffer.append({
                'timestamp_iso': ts_iso,
                'sample_id': sample_id,
                'job_id': job.job_id,
                'job_type': job.job_type,
                'queue': queue_name or '',
                'filepath': os.path.basename(fp),
                'duration_seconds': dur_str,
                'status': status,
                'error': err_short,
                'work_dir': work_dir or ''
            })
            
            # Flush buffer if it reaches threshold (always async to avoid blocking)
            if len(self._csv_buffer) >= self._csv_buffer_size:
                try:
                    asyncio.create_task(self._flush_csv_buffer())
                except Exception:
                    pass
            # If buffer is dangerously large, log warning but don't block
            elif len(self._csv_buffer) >= self._csv_buffer_max_size:
                try:
                    # Create urgent flush task but don't await it
                    asyncio.create_task(self._flush_csv_buffer())
                except Exception:
                    pass
        except Exception:
            pass
        
        # Invalidate stats cache since we're updating stats
        self._stats_cache = None

        # Release backpressure for this job BEFORE awaiting submit_jobs() for downstream
        # triggers. Otherwise, at max_total_waiting (~20480), _dispatch_ready_job blocks
        # forever in _wait_for_global_capacity() and never returns — inflight is not
        # decremented and _drain_waiting_jobs never runs (Coordinator deadlock; GUI frozen).
        try:
            jt_bp = job.job_type
            if self.inflight_by_type.get(jt_bp, 0) > 0:
                self.inflight_by_type[jt_bp] -= 1
                if self.inflight_by_type[jt_bp] == 0:
                    self.inflight_by_type.pop(jt_bp, None)
            qn_bp = queue_name or job_queue_of(jt_bp)
            if self.inflight_by_queue.get(qn_bp, 0) > 0:
                self.inflight_by_queue[qn_bp] -= 1
                if self.inflight_by_queue[qn_bp] == 0:
                    self.inflight_by_queue.pop(qn_bp, None)
            await self._drain_waiting_jobs(qn_bp, jt_bp)
        except Exception:
            pass

        if ok:
            if not is_skipped:
                self.completed_count += 1
                self.completed_by_type[job.job_type] = (
                    self.completed_by_type.get(job.job_type, 0) + 1
                )
            else:
                # Track skip for stats
                self.total_skipped += 1
            # trigger downstream (preferred). For bed_conversion, schedule
            # classifiers/slow once per sample instead of per-file.
            requested = job.context.metadata.get("original_workflow", [])
            req_types = [s.split(":", 1)[1] for s in requested] if requested else None
            if job.job_type == "bed_conversion":
                sample_id = (
                    ctx.get_sample_id() if hasattr(ctx, "get_sample_id") else "unknown"
                )
                # For each requested classification/slow type, allow at most
                # one submitted at a time globally; keep at most one waiting.
                for t in TRIGGERS.get("bed_conversion", []):
                    if req_types is not None and t not in req_types:
                        continue
                    if t not in CLASSIFICATION_TYPES:
                        continue
                    busy = self.classif_pending_by_type.get(t, 0) >= 1
                    if not busy:
                        q = job_queue_of(t)
                        j = Job(next(_job_id_counter), t, q, [f"{q}:{t}"], 0, ctx)
                        # submit directly and mark pending
                        proc = self.processors.get(t)
                        if proc is not None:
                            if getattr(self, "using_pools", False):
                                try:
                                    proc.enqueue.remote(j)
                                except Exception:
                                    pass
                            else:
                                ref = proc.process.remote(j)
                                self._inflight[ref] = j
                            self.classif_pending_by_type[t] = (
                                self.classif_pending_by_type.get(t, 0) + 1
                            )
                            # Record submission for totals used by GUI
                            try:
                                self.submitted_by_type[j.job_type] = (
                                    self.submitted_by_type.get(j.job_type, 0) + 1
                                )
                            except Exception:
                                pass
                            self.total_enqueued += 1
                            # Update per-sample aggregate immediately for GUI samples view
                            try:
                                sid_local = (
                                    ctx.get_sample_id()
                                    if hasattr(ctx, "get_sample_id")
                                    else "unknown"
                                )
                                if sid_local != "unknown":
                                    ent_local = self.samples_by_id.get(sid_local)
                                    if ent_local is None:
                                        ent_local = {
                                            "sample_id": sid_local,
                                            "active_jobs": 0,
                                            "pending_jobs": 0,
                                            "total_jobs": 0,
                                            "completed_jobs": 0,
                                            "failed_jobs": 0,
                                            "job_types": set(),
                                            "last_seen": time.time(),
                                        }
                                        self.samples_by_id[sid_local] = ent_local
                                    ent_local["total_jobs"] += 1
                                    ent_local["active_jobs"] += 1
                                    try:
                                        ent_local["job_types"].add(j.job_type)
                                    except Exception:
                                        pass
                                    ent_local["last_seen"] = time.time()
                            except Exception:
                                pass
                            self.active[j.job_id] = {
                                "job_type": j.job_type,
                                "filepath": j.context.filepath,
                                "queue": q,
                                "start_time": time.time(),
                            }
                    else:
                        # store single waiting job if none recorded yet
                        if self.classif_waiting_by_type.get(t) is None:
                            q = job_queue_of(t)
                            self.classif_waiting_by_type[t] = Job(
                                next(_job_id_counter), t, q, [f"{q}:{t}"], 0, ctx
                            )
            else:
                triggered_jobs: List[Job] = []
                for t in TRIGGERS.get(job.job_type, []):
                    if (req_types is None) or (t in req_types):
                        q = job_queue_of(t)
                        
                        # Only create CNV jobs if preprocessing completed successfully
                        if job.job_type == "preprocessing":
                            if 'preprocessing' not in ctx.results:
                                continue
                            
                            prep_result = ctx.results['preprocessing']
                            prep_status = prep_result.get('status', 'unknown')
                            
                            if prep_status != 'success':
                                continue
                        
                        # Use the updated context from the worker (which has the extracted sample ID) instead of the original context
                        cnv_job = Job(next(_job_id_counter), t, q, [f"{q}:{t}"], 0, ctx)
                        
                        triggered_jobs.append(cnv_job)
                if (not is_skipped) and triggered_jobs:
                    # For CNV jobs, use batching if enabled
                    if self.enable_batching and self.sample_batcher:
                        # Check if any of the triggered jobs are CNV jobs
                        cnv_jobs = [j for j in triggered_jobs if j.job_type == "cnv"]
                        other_jobs = [j for j in triggered_jobs if j.job_type != "cnv"]
                        
                        # Submit non-CNV jobs immediately
                        if other_jobs:
                            await self.submit_jobs(other_jobs)
                        
                        # For CNV jobs, add them to the batcher
                        for cnv_job in cnv_jobs:
                            try:
                                sid_pending = cnv_job.context.get_sample_id() if cnv_job.context else "unknown"
                            except Exception:
                                sid_pending = "unknown"
                            try:
                                if sid_pending != "unknown":
                                    ent_pending = self.samples_by_id.get(sid_pending)
                                    if ent_pending is None:
                                        ent_pending = {
                                            "sample_id": sid_pending,
                                            "active_jobs": 0,
                                            "pending_jobs": 0,
                                            "total_jobs": 0,
                                            "completed_jobs": 0,
                                            "failed_jobs": 0,
                                            "job_types": set(),
                                            "last_seen": time.time(),
                                        }
                                        self.samples_by_id[sid_pending] = ent_pending
                                    ent_pending["pending_jobs"] = ent_pending.get("pending_jobs", 0) + 1
                                    ent_pending["last_seen"] = time.time()
                            except Exception:
                                pass
                            batches = await self.sample_batcher.add_job(cnv_job)
                            if batches:
                                try:
                                    for batch in batches:
                                        sid_batch = batch.sample_id
                                        if sid_batch != "unknown":
                                            ent_pending = self.samples_by_id.get(sid_batch)
                                            if ent_pending and ent_pending.get("pending_jobs", 0) > 0:
                                                ent_pending["pending_jobs"] = max(
                                                    0,
                                                    ent_pending.get("pending_jobs", 0) - len(batch.contexts),
                                                )
                                                ent_pending["last_seen"] = time.time()
                                except Exception:
                                    pass
                                # Convert BatchedJob to regular Job for processing
                                regular_jobs = self._convert_batched_jobs(batches)
                                await self._submit_jobs_internal(regular_jobs)
                    else:
                        # No batching, submit all jobs normally
                        await self.submit_jobs(triggered_jobs)
                else:
                    # advance linear chain if any
                    nxt = job.next_job()
                    if (not is_skipped) and nxt:
                        await self.submit_jobs([nxt])
        else:
            # Only record as failed if not a deliberate skip
            if not is_skipped:
                self.failed_count += 1
                self.failed_by_type[job.job_type] = (
                    self.failed_by_type.get(job.job_type, 0) + 1
                )
                # Record failed job for tracking and OOM diagnosis
                try:
                    err_str = (err or "")[:1500]  # Keep more of traceback so root cause is visible
                    err_lower = err_str.lower()
                    # Ray OOM: worker is killed by raylet; error often says "worker died" / "connection" / "SIGKILL" not "oom"
                    is_oom = (
                        "oom" in err_lower
                        or "out of memory" in err_lower
                        or "outofmemory" in err_lower
                        or "worker killed" in err_lower
                        or "killed (oom)" in err_lower
                        or "memory error" in err_lower
                        or "worker died" in err_lower
                        or "worker unexpectedly exits" in err_lower
                        or "worker exit" in err_lower
                        or "sigkill" in err_lower
                        or "connection error" in err_lower
                        or "system_error" in err_lower
                    )
                    sid_fail = (
                        job.context.get_sample_id()
                        if hasattr(job.context, "get_sample_id")
                        else "unknown"
                    )
                    fp_fail = getattr(job.context, "filepath", "") or ""
                    self._failed_jobs_log.append({
                        "sample_id": sid_fail,
                        "filepath": fp_fail,
                        "job_type": job.job_type,
                        "error": err_str,
                        "is_oom": is_oom,
                        "failure_kind": failure_kind or "",
                    })
                except Exception:
                    pass

        # per-sample finish updates
        try:
            sid2 = (
                job.context.get_sample_id()
                if hasattr(job.context, "get_sample_id")
                else "unknown"
            )
            ent2 = self.samples_by_id.get(sid2)
            if ent2 is not None:
                if ent2.get("active_jobs", 0) > 0:
                    ent2["active_jobs"] -= 1
                # Do not count deliberately skipped preprocessing as completed or failed
                is_skipped = False
                try:
                    is_skipped = (job.job_type == "preprocessing") and bool(
                        job.context.metadata.get("skip_downstream", False)
                    )
                except Exception:
                    is_skipped = False
                if not is_skipped:
                    if ok:
                        ent2["completed_jobs"] = ent2.get("completed_jobs", 0) + 1
                    else:
                        ent2["failed_jobs"] = ent2.get("failed_jobs", 0) + 1
                ent2["last_seen"] = time.time()
        except Exception:
            pass
        
        # Periodically clean up completed samples to prevent memory growth
        self._jobs_since_last_cleanup += 1
        if self._jobs_since_last_cleanup >= self._sample_cleanup_interval:
            self._jobs_since_last_cleanup = 0
            try:
                asyncio.create_task(self._cleanup_completed_samples())
            except Exception:
                pass
        
        # Periodically clean up CNV cache to prevent memory growth
        self._jobs_since_cnv_cleanup += 1
        if self._jobs_since_cnv_cleanup >= self._cnv_cache_cleanup_interval:
            self._jobs_since_cnv_cleanup = 0
            try:
                asyncio.create_task(self._cleanup_cnv_cache())
            except Exception:
                pass

        # Promote next waiting job for this (type, sample) if any (CNV-only if enabled)
        if job.job_type in SERIALIZE_BY_TYPE_PER_SAMPLE:
            sid = job.context.get_sample_id() or "unknown"
            key_ts = (job.job_type, sid)
            if key_ts in self.running_by_type_sample:
                self.running_by_type_sample[key_ts] = max(
                    0, self.running_by_type_sample.get(key_ts, 0) - 1
                )
                if self.running_by_type_sample[key_ts] == 0:
                    self.running_by_type_sample.pop(key_ts, None)
            queue_list = self.waiting_by_type_sample.get(key_ts, [])
            if queue_list:
                nxt_job = queue_list.pop(0)
                if not queue_list:
                    self.waiting_by_type_sample.pop(key_ts, None)
                self.pending_by_type_sample[key_ts] = max(
                    0, self.pending_by_type_sample.get(key_ts, 0) - 1
                )
                if self.pending_by_type_sample[key_ts] == 0:
                    self.pending_by_type_sample.pop(key_ts, None)
                try:
                    if sid != "unknown":
                        ent_pending = self.samples_by_id.get(sid)
                        if ent_pending and ent_pending.get("pending_jobs", 0) > 0:
                            ent_pending["pending_jobs"] -= 1
                            ent_pending["last_seen"] = time.time()
                except Exception:
                    pass
                self.running_by_type_sample[key_ts] = (
                    self.running_by_type_sample.get(key_ts, 0) + 1
                )
                proc2 = self.processors.get(nxt_job.job_type)
                if proc2 is not None:
                    ref2 = proc2.process.remote(nxt_job)
                    self._inflight[ref2] = nxt_job
                    # Record submission for totals used by GUI
                    try:
                        self.submitted_by_type[nxt_job.job_type] = (
                            self.submitted_by_type.get(nxt_job.job_type, 0) + 1
                        )
                    except Exception:
                        pass
                    self.total_enqueued += 1
                    # Update per-sample aggregate immediately for GUI samples view
                    try:
                        sid_local2 = (
                            nxt_job.context.get_sample_id()
                            if hasattr(nxt_job.context, "get_sample_id")
                            else "unknown"
                        )
                        if sid_local2 != "unknown":
                            ent_local2 = self.samples_by_id.get(sid_local2)
                            if ent_local2 is None:
                                ent_local2 = {
                                    "sample_id": sid_local2,
                                    "active_jobs": 0,
                                    "pending_jobs": 0,
                                    "total_jobs": 0,
                                    "completed_jobs": 0,
                                    "failed_jobs": 0,
                                    "job_types": set(),
                                    "last_seen": time.time(),
                                }
                                self.samples_by_id[sid_local2] = ent_local2
                            ent_local2["total_jobs"] += 1
                            ent_local2["active_jobs"] += 1
                            try:
                                ent_local2["job_types"].add(nxt_job.job_type)
                            except Exception:
                                pass
                            ent_local2["last_seen"] = time.time()
                    except Exception:
                        pass
                    self.active[nxt_job.job_id] = {
                        "job_type": nxt_job.job_type,
                        "filepath": nxt_job.context.filepath,
                        "queue": job_queue_of(nxt_job.job_type),
                        "start_time": time.time(),
                    }

        if job.job_type in CLASSIFICATION_TYPES:
            # decrement pending for this type
            if self.classif_pending_by_type.get(job.job_type, 0) > 0:
                self.classif_pending_by_type[job.job_type] -= 1
                if self.classif_pending_by_type[job.job_type] == 0:
                    self.classif_pending_by_type.pop(job.job_type, None)
            nxt = self.classif_waiting_by_type.get(job.job_type)
            if nxt is not None:
                self.classif_waiting_by_type[job.job_type] = None
                proc3 = self.processors.get(job.job_type)
                if proc3 is not None:
                    if getattr(self, "using_pools", False):
                        try:
                            proc3.enqueue.remote(nxt)
                        except Exception:
                            pass
                    else:
                        ref3 = proc3.process.remote(nxt)
                        self._inflight[ref3] = nxt
                    # Record submission for totals used by GUI
                    try:
                        self.submitted_by_type[job.job_type] = (
                            self.submitted_by_type.get(job.job_type, 0) + 1
                        )
                    except Exception:
                        pass
                    self.classif_pending_by_type[job.job_type] = 1
                    self.total_enqueued += 1
                    # Update per-sample aggregate immediately for GUI samples view
                    try:
                        sid_local3 = (
                            nxt.context.get_sample_id()
                            if hasattr(nxt, "context")
                            and hasattr(nxt.context, "get_sample_id")
                            else "unknown"
                        )
                        if sid_local3 != "unknown":
                            ent_local3 = self.samples_by_id.get(sid_local3)
                            if ent_local3 is None:
                                ent_local3 = {
                                    "sample_id": sid_local3,
                                    "active_jobs": 0,
                                    "pending_jobs": 0,
                                    "total_jobs": 0,
                                    "completed_jobs": 0,
                                    "failed_jobs": 0,
                                    "job_types": set(),
                                    "last_seen": time.time(),
                                }
                                self.samples_by_id[sid_local3] = ent_local3
                            ent_local3["total_jobs"] += 1
                            ent_local3["active_jobs"] += 1
                            try:
                                ent_local3["job_types"].add(nxt.job_type)
                            except Exception:
                                pass
                            ent_local3["last_seen"] = time.time()
                    except Exception:
                        pass
                    self.active[nxt.job_id] = {
                        "job_type": nxt.job_type,
                        "filepath": nxt.context.filepath,
                        "queue": job_queue_of(nxt.job_type),
                        "start_time": time.time(),
                    }

    async def _drain_loop(self):
        while True:
            try:
                if not self._inflight:
                    await asyncio.sleep(0.1)
                    continue
                now = time.time()
                timeout_sec = getattr(self, "job_timeout_seconds", 0) or 0
                # Cancel any jobs that have exceeded the per-job timeout
                if timeout_sec > 0:
                    for ref, job in list(self._inflight.items()):
                        start = None
                        try:
                            info = self.active.get(job.job_id)
                            if isinstance(info, dict):
                                start = info.get("start_time")
                        except Exception:
                            pass
                        if start is not None and (now - start) >= timeout_sec:
                            try:
                                ray.cancel(ref)
                            except Exception:
                                pass
                            self._inflight.pop(ref, None)
                            err_msg = f"Job timed out after {timeout_sec}s"
                            await self._on_finish(job, False, job.context, err_msg, "timed_out")
                refs = list(self._inflight.keys())
                if not refs:
                    await asyncio.sleep(0.1)
                    continue
                ready, _ = ray.wait(refs, num_returns=1, timeout=0.1)
                for r in ready:
                    job = self._inflight.pop(r, None)
                    ok, ctx, err = True, None, None
                    try:
                        ctx = await r
                    except Exception as e:
                        ok, err = False, str(e)
                        ctx = job.context if job else WorkflowContext(filepath="")
                    if job:
                        await self._on_finish(job, ok, ctx, err)
            except Exception:
                await asyncio.sleep(0.1)
    
    async def _batch_timeout_loop(self):
        """Periodically check for timed-out batches"""
        while not self._shutdown_requested:
            try:
                # Check for timed-out batches
                timed_out_batches = await self.sample_batcher.check_timeouts()
                
                # Submit timed-out batches
                if timed_out_batches:
                    try:
                        for batch in timed_out_batches:
                            sid_batch = batch.sample_id
                            if sid_batch != "unknown":
                                ent_pending = self.samples_by_id.get(sid_batch)
                                if ent_pending and ent_pending.get("pending_jobs", 0) > 0:
                                    ent_pending["pending_jobs"] = max(
                                        0,
                                        ent_pending.get("pending_jobs", 0) - len(batch.contexts),
                                    )
                                    ent_pending["last_seen"] = time.time()
                    except Exception:
                        pass
                    regular_jobs = self._convert_batched_jobs(timed_out_batches)
                    await self._submit_jobs_internal(regular_jobs)
                
                await asyncio.sleep(1.0)  # Check every second
                
            except Exception:
                await asyncio.sleep(1.0)
    
    def _convert_batched_jobs(self, batched_jobs: List[BatchedJob]) -> List[Job]:
        """Convert BatchedJob objects to regular Job objects for processing"""
        regular_jobs = []
        for batched_job in batched_jobs:
            # Create a regular job that will handle the batch
            regular_job = Job(
                job_id=batched_job.job_id,
                job_type=batched_job.job_type,
                origin=batched_job.origin,
                workflow=batched_job.workflow,
                step=batched_job.step,
                context=batched_job.contexts[0]  # Use first context as primary
            )
            # Store batch information in metadata
            regular_job.context.metadata["_batched_job"] = batched_job
            regular_jobs.append(regular_job)
        
        return regular_jobs
    
    def register_batched_handler(self, job_type: str, handler: Callable[[BatchedJob], None]) -> None:
        """Register a handler that can process batched jobs"""
        # Wrap the handler to extract BatchedJob from regular Job
        def wrapped_handler(job: Job) -> None:
            batched_job = job.context.metadata.get("_batched_job")
            if batched_job:
                handler(batched_job)
            else:
                # Fallback to single file processing
                single_context = job.context
                single_batch = BatchedJob(
                    job_id=job.job_id,
                    job_type=job.job_type,
                    origin=job.origin,
                    workflow=job.workflow,
                    step=job.step,
                    contexts=[single_context],
                    batch_id=f"single_{job.job_id}",
                    sample_id=single_context.get_sample_id()
                )
                handler(single_batch)
        
        # Register the wrapped handler
        self.register_handler(job_type, wrapped_handler)

    async def mark_running(self, job: Job):
        if job.job_type in DEDUP_TYPES:
            key = (job.context.get_sample_id(), job.job_type)
            self.pending[key] = max(0, self.pending.get(key, 0) - 1)
            if self.pending[key] == 0:
                self.pending.pop(key, None)
            self.running[key] = self.running.get(key, 0) + 1

    async def mark_handler_started(self, job_id: int) -> None:
        """Worker calls this when the real handler body begins (Ray task running)."""
        try:
            info = self.active.get(int(job_id))
            if info is None:
                return
            if info.get("handler_start_time") is not None:
                return
            info["handler_start_time"] = time.time()
        except Exception:
            pass

    async def stats(self) -> Dict[str, Any]:
        # Check cache first for performance
        now = time.time()
        if self._stats_cache and (now - self._stats_cache_time) < self._stats_cache_ttl:
            return self._stats_cache
        
        # per-queue active summary (also expose as 'active_by_worker' for GUI compat)
        # IMPORTANT: We keep two representations:
        # - active_by_queue: truncated to the 2 most recently started jobs (UI display)
        # - active_count_by_queue/job_type: full counts (debugging "Act" discrepancies)
        active_by_queue: Dict[str, List[Dict[str, Any]]] = {}
        active_count_by_queue: Dict[str, int] = {}
        active_count_by_job_type: Dict[str, int] = {}
        for jid, info in self.active.items():
            q = info["queue"]
            jt = info["job_type"]
            active_count_by_queue[q] = active_count_by_queue.get(q, 0) + 1
            active_count_by_job_type[jt] = active_count_by_job_type.get(jt, 0) + 1
            hs = info.get("handler_start_time")
            st = info.get("start_time")
            if hs is not None:
                duration_s = int(now - float(hs))
                dispatch_wait = (
                    int(float(hs) - float(st)) if st is not None else None
                )
            else:
                duration_s = int(now - float(st)) if st is not None else 0
                dispatch_wait = None
            active_by_queue.setdefault(q, []).append(
                {
                    "job_id": jid,
                    "job_type": jt,
                    "filepath": info["filepath"],
                    "duration": duration_s,
                    "dispatch_wait_seconds": dispatch_wait,
                    "start_time": info["start_time"],
                }
            )

        # Sort each queue's list by start_time descending (newest first), then keep only first 2
        for q in active_by_queue:
            active_by_queue[q] = sorted(
                active_by_queue[q], key=lambda x: x["start_time"], reverse=True
            )[:2]
            for entry in active_by_queue[q]:
                entry.pop("start_time", None)

        # Fallback path: if coordinator active bookkeeping is empty/missing while
        # queue actors still report running jobs, surface those in monitor stats.
        try:
            runtime_rows = await asyncio.gather(
                *[proc.runtime_status.remote() for proc in self.processors.values()],
                return_exceptions=True,
            )
        except Exception:
            runtime_rows = []
        runtime_active_count_by_queue: Dict[str, int] = {}
        runtime_active_count_by_job_type: Dict[str, int] = {}
        runtime_active_by_queue: Dict[str, List[Dict[str, Any]]] = {}
        for row in runtime_rows:
            if isinstance(row, Exception) or not isinstance(row, dict):
                continue
            qn = str(row.get("queue") or "")
            rc = int(row.get("running_count") or 0)
            if qn:
                runtime_active_count_by_queue[qn] = runtime_active_count_by_queue.get(qn, 0) + rc
            for j in (row.get("inflight") or []):
                if not isinstance(j, dict):
                    continue
                jt = str(j.get("job_type") or "")
                if jt:
                    runtime_active_count_by_job_type[jt] = runtime_active_count_by_job_type.get(jt, 0) + 1
                if qn:
                    runtime_active_by_queue.setdefault(qn, []).append(
                        {
                            "job_id": j.get("job_id"),
                            "job_type": jt,
                            "filepath": j.get("filepath", ""),
                            "duration": int(j.get("duration", 0) or 0),
                            "dispatch_wait_seconds": None,
                        }
                    )

        if len(self.active) == 0:
            # Use runtime fallback only when primary active map is empty.
            if runtime_active_count_by_queue:
                active_count_by_queue = runtime_active_count_by_queue
            if runtime_active_count_by_job_type:
                active_count_by_job_type = runtime_active_count_by_job_type
            if runtime_active_by_queue:
                active_by_queue = {
                    q: sorted(v, key=lambda x: x.get("duration", 0), reverse=True)[:2]
                    for q, v in runtime_active_by_queue.items()
                }

        # A compact dump of the longest-running active jobs (useful for "stuck" debugging)
        active_top_by_duration = sorted(
            [
                {
                    "job_id": jid,
                    "job_type": info["job_type"],
                    "queue": info["queue"],
                    "filepath": info["filepath"],
                    "duration": (
                        int(now - float(info["handler_start_time"]))
                        if info.get("handler_start_time") is not None
                        else int(now - float(info["start_time"]))
                    ),
                }
                for jid, info in self.active.items()
            ],
            key=lambda x: x["duration"],
            reverse=True,
        )[:10]

        # build category counts
        def _cat_of(jt: str) -> str:
            if jt == "preprocessing":
                return "preprocessing"
            if jt in {"mgmt", "cnv", "target", "fusion"}:
                return jt
            if jt in CLASSIFICATION_TYPES:
                return "classification"
            return "other"

        running_counts: Dict[str, int] = {
            "preprocessing": 0,
            "mgmt": 0,
            "cnv": 0,
            "target": 0,
            "fusion": 0,
            "classification": 0,
            "other": 0,
        }
        totals_counts: Dict[str, int] = {
            "preprocessing": 0,
            "mgmt": 0,
            "cnv": 0,
            "target": 0,
            "fusion": 0,
            "classification": 0,
            "other": 0,
        }
        for info in self.active.values():
            running_counts[_cat_of(info["job_type"])] += 1
        # totals are the submitted (ever) per type/category
        for jt, c in self.submitted_by_type.items():
            totals_counts[_cat_of(jt)] += int(c or 0)
        # Count waiting (serialized) jobs so the monitor can account for them
        try:
            waiting_serialized = sum(
                len(v) for v in getattr(self, "waiting_by_type_sample", {}).values()
            )
        except Exception:
            waiting_serialized = 0
        # Count global backpressure queues (per-type, per-queue, and global)
        try:
            waiting_global = sum(
                len(v) for v in getattr(self, "waiting_by_type_global", {}).values()
            )
        except Exception:
            waiting_global = 0
        try:
            waiting_global += sum(
                len(v) for v in getattr(self, "waiting_by_queue", {}).values()
            )
        except Exception:
            pass
        try:
            waiting_global += len(getattr(self, "waiting_global", []))
        except Exception:
            pass

        waiting_count_by_queue: Dict[str, int] = {}
        try:
            # For debug: how much is queued per queue name
            for qn, lst in getattr(self, "waiting_by_queue", {}).items():
                waiting_count_by_queue[qn] = len(lst)
        except Exception:
            waiting_count_by_queue = {}

        waiting_count_by_type_global: Dict[str, int] = {}
        try:
            # For debug: how much is queued in the per-job-type global backpressure lists
            for jt, lst in getattr(self, "waiting_by_type_global", {}).items():
                waiting_count_by_type_global[jt] = len(lst)
        except Exception:
            waiting_count_by_type_global = {}

        inflight_by_type_debug: Dict[str, int] = {}
        inflight_by_queue_debug: Dict[str, int] = {}
        try:
            inflight_by_type_debug = dict(getattr(self, "inflight_by_type", {}) or {})
        except Exception:
            inflight_by_type_debug = {}
        try:
            inflight_by_queue_debug = dict(getattr(self, "inflight_by_queue", {}) or {})
        except Exception:
            inflight_by_queue_debug = {}
        # Samples payload for GUI
        samples_payload: List[Dict[str, Any]] = []
        try:
            for sid, ent in self.samples_by_id.items():
                sample_data = {
                    "sample_id": sid,
                    "active_jobs": ent.get("active_jobs", 0),
                    "pending_jobs": ent.get("pending_jobs", 0),
                    "total_jobs": ent.get("total_jobs", 0),
                    "completed_jobs": ent.get("completed_jobs", 0),
                    "failed_jobs": ent.get("failed_jobs", 0),
                    "job_types": list(ent.get("job_types", set())),
                    "last_seen": ent.get("last_seen", now),
                }
                samples_payload.append(sample_data)
        except Exception as e:
            logging.error(f"[Ray] Error building samples payload: {e}")
            samples_payload = []

        result = {
            "completed": self.completed_count,
            "failed": self.failed_count,
            "total_enqueued": self.total_enqueued,
            "total_skipped": self.total_skipped,
            "completed_by_type": dict(self.completed_by_type),
            "failed_by_type": dict(self.failed_by_type),
            "active_by_queue": active_by_queue,
            "active_by_worker": active_by_queue,  # compatibility for GUI table
            "active_count": len(self.active),
            "active_count_by_queue": active_count_by_queue,
            "active_count_by_job_type": active_count_by_job_type,
            "active_top_by_duration": active_top_by_duration,
            "waiting_serialized": waiting_serialized,
            "waiting_global": waiting_global,
            "waiting_count_by_queue": waiting_count_by_queue,
            "waiting_count_by_type_global": waiting_count_by_type_global,
            "inflight_by_type_debug": inflight_by_type_debug,
            "inflight_by_queue_debug": inflight_by_queue_debug,
            "running_by_category": running_counts,
            "totals_by_category": totals_counts,
            "samples": samples_payload,
            "runtime_active_count_by_queue": runtime_active_count_by_queue,
            "runtime_active_count_by_job_type": runtime_active_count_by_job_type,
            # Failed job tracking: last 500 entries for inspection / OOM diagnosis
            "failed_jobs_list": list(self._failed_jobs_log)[-500:],
            "oom_count": sum(1 for e in self._failed_jobs_log if e.get("is_oom")),
        }
        if result["active_count"] == 0 and runtime_active_count_by_queue:
            result["active_count"] = int(sum(runtime_active_count_by_queue.values()))
        
        # Cache the result for future calls
        self._stats_cache = result
        self._stats_cache_time = now
        
        return result

    async def get_cnv_cache_stats(self) -> Dict[str, Any]:
        """
        Get CNV cache statistics for monitoring purposes.
        
        Returns:
            Dictionary with CNV cache statistics
        """
        if get_sample_cache_stats is None:
            return {"error": "CNV module not available"}
        
        try:
            return get_sample_cache_stats()
        except Exception as e:
            return {"error": str(e)}

    async def shutdown(self) -> bool:
        """
        Gracefully shutdown the coordinator and all its actors.
        
        Returns:
            True if shutdown was successful
        """
        try:
            print("[SHUTDOWN] Stopping workflow coordinator...")
            # Stop accepting new jobs
            self._shutdown_requested = True
            print("[SHUTDOWN] Stopped accepting new jobs")

            # Write failed jobs log so user can see which jobs failed (e.g. OOM)
            try:
                path = self._write_failed_jobs_file()
                if path:
                    print(f"[SHUTDOWN] Wrote failed jobs log: {path}")
            except Exception as e:
                print(f"[SHUTDOWN] Warning: Could not write failed jobs log: {e}")
            
            # Finalize target accumulation for all samples with batch BAMs pending merge
            try:
                print("[SHUTDOWN] Finalizing target.bam files for all samples...")
                await self._finalize_target_bams()
                print("[SHUTDOWN] Target BAM finalization complete")
            except Exception as e:
                print(f"[SHUTDOWN] Warning: Error during target BAM finalization: {e}")
            
            # Flush any remaining CSV buffer entries
            try:
                print("[SHUTDOWN] Flushing CSV buffer entries...")
                await self._flush_csv_buffer()
                print("[SHUTDOWN] CSV buffer flushed")
            except Exception as e:
                print(f"[SHUTDOWN] Warning: Error flushing CSV buffer: {e}")
            
            # Cancel all inflight tasks
            if hasattr(self, '_inflight') and self._inflight:
                inflight_count = len(self._inflight)
                if inflight_count > 0:
                    print(f"[SHUTDOWN] Cancelling {inflight_count} inflight tasks...")
                    for ref in list(self._inflight.keys()):
                        try:
                            ray.cancel(ref)
                        except Exception:
                            pass
                    print("[SHUTDOWN] Inflight tasks cancelled")
            
            # Shutdown all pool actors gracefully first
            if hasattr(self, 'queues') and self.queues:
                total_actors = sum(len(actors) for actors in self.queues.values())
                if total_actors > 0:
                    print(f"[SHUTDOWN] Shutting down {total_actors} pool actors...")
                    for queue_name, actors in self.queues.items():
                        for actor in actors:
                            try:
                                await actor.shutdown.remote()
                            except Exception:
                                pass
                            try:
                                ray.kill(actor)
                            except Exception:
                                pass
                    print("[SHUTDOWN] Pool actors shut down")
            
            # Kill all processor actors
            if hasattr(self, 'processors') and self.processors:
                processor_count = len(self.processors)
                if processor_count > 0:
                    print(f"[SHUTDOWN] Shutting down {processor_count} processor actors...")
                    for processor in self.processors.values():
                        try:
                            ray.kill(processor)
                        except Exception:
                            pass
                    print("[SHUTDOWN] Processor actors shut down")
            
            print("[SHUTDOWN] Coordinator shutdown complete")
            return True
        except Exception as e:
            print(f"[SHUTDOWN] Error during shutdown: {e}")
            return False


@ray.remote
class Pool:
    def __init__(
        self,
        queue_name: str,
        max_parallel: int = 1,
        job_timeout_seconds: int = 0,
    ):
        self.queue_name = queue_name
        self.max_parallel = max(1, int(max_parallel))
        self.job_timeout_seconds = max(0, int(job_timeout_seconds))
        # job_type -> (remote_func, resource_options)
        self.handlers: Dict[str, Tuple[Any, Dict[str, Any]]] = {}
        # callback to Coordinator actor
        self._coordinator = None
        self._coordinator_name: Optional[str] = None
        # ref -> (job, start_time) for timeout and completion
        self._inflight_data: Dict[Any, Tuple[Job, float]] = {}
        self._timed_out_jobs: Set[str] = set()
        self._running_count: int = 0
        self._shutdown_requested: bool = False
        self._timeout_task: Optional[asyncio.Task] = None
        
        # Initialize memory manager for long-running pool actor
        self.memory_manager = None
        if MemoryManager is not None and not DISABLE_MEMORY_MANAGER:
            try:
                # Configure memory management based on queue type
                gc_every = 30 if queue_name in {"analysis", "classif"} else 50
                rss_trigger = 1536 if queue_name in {"analysis", "classif"} else 2048
                restart_every = 7500 if queue_name in {"analysis", "classif"} else 10000
                restart_rss_trigger = 3072 if queue_name in {"analysis", "classif"} else 4096
                self.memory_manager = MemoryManager(
                    gc_every=gc_every,
                    rss_trigger_mb=rss_trigger,
                    enable_malloc_trim=True,
                    restart_every=restart_every,
                    restart_rss_trigger_mb=restart_rss_trigger
                )
            except Exception as e:
                # Fallback if memory manager fails to initialize
                pass

    async def set_coordinator_name(self, coord_name: str):
        # store Coordinator actor name; resolve handle lazily
        self._coordinator_name = coord_name
        try:
            self._coordinator = ray.get_actor(coord_name)
        except Exception:
            self._coordinator = None
        return True

    def _get_coordinator(self):
        if self._coordinator is None and self._coordinator_name:
            try:
                self._coordinator = ray.get_actor(self._coordinator_name)
            except Exception:
                self._coordinator = None
        return self._coordinator
    
    def should_restart(self) -> bool:
        """Check if actor should restart due to memory management."""
        if self.memory_manager is not None:
            return self.memory_manager.is_restart_requested()
        return False
    
    def get_restart_reason(self) -> Optional[str]:
        """Get the reason for restart request."""
        if self.memory_manager is not None:
            return self.memory_manager.get_restart_reason()
        return None

    async def register_handler(
        self, job_type: str, remote_func, resource_options: Dict[str, Any]
    ):
        self.handlers[job_type] = (remote_func, resource_options)

    async def runtime_status(self) -> Dict[str, Any]:
        """
        Lightweight runtime status for progress fallback/debugging.
        Returns queue-local running count and a short list of inflight jobs.
        """
        now = time.time()
        inflight: List[Dict[str, Any]] = []
        try:
            for ref, (job, start_time) in list(self._inflight_data.items())[:4]:
                try:
                    inflight.append(
                        {
                            "job_id": job.job_id,
                            "job_type": job.job_type,
                            "filepath": job.context.filepath if job.context else "",
                            "duration": int(now - float(start_time)),
                        }
                    )
                except Exception:
                    continue
        except Exception:
            inflight = []
        return {
            "queue": self.queue_name,
            "running_count": int(self._running_count),
            "inflight": inflight,
        }

    async def enqueue(self, job: Job):
        # Check if shutdown has been requested
        if self._shutdown_requested:
            return
        
        # tell coordinator that we're starting this job (updates dedup maps)
        # We cannot call back into Coordinator directly from here without a reference;
        # pass coordinator handle via callback's bound actor method (Ray passes actor ref implicitly).
        # Instead, we expose a small trampoline: the driver calls mark_running before enqueue.
        # For simplicity, Coordinator calls mark_running() before Pool.enqueue().
        tup = self.handlers.get(job.job_type)
        if tup is None:
            ctx = job.context
            ctx.add_error(
                job.job_type,
                f"No handler for {job.job_type} in queue {self.queue_name}",
            )
            coord = self._get_coordinator()
            if coord is not None:
                await coord._on_finish.remote(
                    job, False, ctx, f"no_handler:{job.job_type}"
                )
            return
        remote_func, opts = tup
        # Apply resource hints at call site
        # Bound concurrency to max_parallel to emulate thread runner
        while self._running_count >= self.max_parallel:
            await asyncio.sleep(0.005)
        ref = remote_func.options(**opts).remote(job)
        self._inflight_data[ref] = (job, time.time())
        self._running_count += 1

        # Start timeout loop once if per-job timeout is enabled
        if self.job_timeout_seconds > 0 and (self._timeout_task is None or self._timeout_task.done()):
            try:
                self._timeout_task = asyncio.create_task(self._timeout_loop())
            except Exception:
                pass

        async def _wait(r):
            ok, ctx, err = True, None, None
            failure_kind = None
            try:
                ctx = await r
                # Treat handler "error return" as failure: context has errors for this job_type
                errors = getattr(ctx, "errors", []) or []
                job_errors = [e for e in errors if isinstance(e, dict) and e.get("job_type") == job.job_type]
                if job_errors:
                    ok = False
                    err = job_errors[-1].get("error", "handler reported error")
                    failure_kind = "returned_error"
            except Exception as e:
                ok, err = False, str(e)
                ctx = job.context
                ctx.add_error(job.job_type, err)
                failure_kind = "exception"
            self._inflight_data.pop(r, None)
            coord2 = self._get_coordinator()
            if coord2 is not None and job.job_id not in self._timed_out_jobs:
                await coord2._on_finish.remote(job, ok, ctx, err, failure_kind)
            elif job.job_id in self._timed_out_jobs:
                self._timed_out_jobs.discard(job.job_id)
            if self._running_count > 0:
                self._running_count -= 1
            
            # Trigger memory cleanup after job completion
            if self.memory_manager is not None:
                try:
                    is_idle = self._running_count == 0
                    cleanup_stats = self.memory_manager.check_and_cleanup(is_idle=is_idle)
                    # Log memory cleanup if it was triggered
                    if cleanup_stats.get('gc_triggered', False):
                        pass  # Memory cleanup performed
                    
                    # Actor restart mechanism disabled - no restart checks
                except Exception as e:
                    # Silently continue if memory management fails
                    pass

        asyncio.create_task(_wait(ref))

    async def shutdown(self) -> bool:
        """
        Gracefully shutdown the pool actor.
        
        Returns:
            True if shutdown was successful
        """
        try:
            # Stop accepting new jobs
            self._shutdown_requested = True
            
            # Cancel all inflight tasks
            if self._inflight_data:
                for ref in list(self._inflight_data):
                    try:
                        ray.cancel(ref)
                    except Exception:
                        pass
            
            return True
        except Exception:
            return False

    async def _timeout_loop(self) -> None:
        """Cancel jobs that exceed job_timeout_seconds and report them as failed."""
        while not self._shutdown_requested and self.job_timeout_seconds > 0:
            try:
                await asyncio.sleep(10)
                now = time.time()
                coord = self._get_coordinator()
                if coord is None:
                    continue
                for ref, (job, start_time) in list(self._inflight_data.items()):
                    if (now - start_time) >= self.job_timeout_seconds:
                        self._timed_out_jobs.add(job.job_id)
                        self._inflight_data.pop(ref, None)
                        if self._running_count > 0:
                            self._running_count -= 1
                        try:
                            ray.cancel(ref)
                        except Exception:
                            pass
                        err_msg = (
                            f"Job timed out after {self.job_timeout_seconds}s"
                        )
                        await coord._on_finish.remote(
                            job, False, job.context, err_msg, "timed_out"
                        )
            except asyncio.CancelledError:
                break
            except Exception:
                pass

    # Adapter so callers using TypeProcessor-style .process() keep working
    async def process(self, job: Job):
        return await self.enqueue(job)


# ---------- Classifier & Runner ----------


def default_file_classifier(filepath: str, plan: List[str], target_panel: str) -> List[Job]:
    ctx = WorkflowContext(filepath)
    ctx.add_metadata("filename", os.path.basename(filepath))
    ctx.add_metadata("created", time.time())
    ctx.add_metadata("target_panel", target_panel)  # Add panel metadata
    try:
        ctx.add_metadata("file_size", os.path.getsize(filepath))
    except OSError:
        ctx.add_metadata("file_size", 0)
    ctx.add_metadata("original_workflow", plan)

    job_id = next(_job_id_counter)
    first = plan[0]
    if ":" not in first:
        raise ValueError(f"Invalid workflow step: {first}")
    q, jt = first.split(":", 1)
    return [Job(job_id, jt, q, plan, 0, ctx)]


def _matches_any_pattern(p: Path, patterns: Optional[List[str]]) -> bool:
    if not patterns:
        return True
    try:
        return any(p.match(pat) for pat in patterns)
    except Exception:
        return True


def _matches_no_ignores(p: Path, ignore_patterns: Optional[List[str]]) -> bool:
    # Always ignore hidden dotfiles/directories by default.
    # This keeps `robin workflow` from getting noisy when watch directories contain
    # files like `.DS_Store` or editor artifacts.
    if p.name.startswith("."):
        return False
    if not ignore_patterns:
        return True
    try:
        return not any(p.match(pat) for pat in ignore_patterns)
    except Exception:
        return True


def _is_bam_path(p: Path) -> bool:
    try:
        return p.suffix.lower() == ".bam"
    except Exception:
        return False


def _bam_filename_kind(filepath: str) -> str:
    """Classify BAM filenames by substring: 'pass' | 'fail' | 'unknown'."""
    try:
        n = os.path.basename(filepath).lower()
    except Exception:
        return "unknown"
    if "pass" in n:
        return "pass"
    if "fail" in n:
        return "fail"
    return "unknown"


def _detect_fail_only_bam_submission(
    paths: List[str],
    patterns: Optional[List[str]],
    ignore_patterns: Optional[List[str]],
    recursive: bool,
) -> bool:
    """
    True iff this submission contains BAMs and *all* BAM filenames contain 'fail'
    and none contain 'pass'. Used to suppress expected downstream failures on
    fail-only read sets.
    """
    bam_total = 0
    bam_fail = 0
    bam_pass = 0
    for p in paths:
        pth = Path(p)
        if pth.is_file():
            files = [pth]
        elif pth.is_dir():
            try:
                files = list(pth.rglob("*") if recursive else pth.glob("*"))
            except Exception:
                files = []
        else:
            files = []
        for f in files:
            try:
                if not f.is_file():
                    continue
                if not _matches_any_pattern(f, patterns) or not _matches_no_ignores(
                    f, ignore_patterns
                ):
                    continue
                if not _is_bam_path(f):
                    continue
                bam_total += 1
                k = _bam_filename_kind(str(f))
                if k == "fail":
                    bam_fail += 1
                elif k == "pass":
                    bam_pass += 1
            except Exception:
                continue
    return bool(bam_total and bam_pass == 0 and bam_fail == bam_total)


async def submit_existing_paths(
    coord,
    paths: List[str],
    plan: List[str],
    patterns: Optional[List[str]] = None,
    ignore_patterns: Optional[List[str]] = None,
    recursive: bool = True,
    work_dir: Optional[str] = None,
) -> None:
    # Fetch reference and target_panel once before processing files (performance optimization)
    coord_reference = None
    coord_target_panel = None
    try:
        coord_reference = await coord.get_reference.remote()
    except Exception:
        pass
    try:
        coord_target_panel = await coord.get_target_panel.remote()
    except Exception:
        pass
    
    fail_only_bam_submission = _detect_fail_only_bam_submission(
        paths=paths,
        patterns=patterns,
        ignore_patterns=ignore_patterns,
        recursive=recursive,
    )

    seed_jobs: List[Job] = []
    for p in paths:
        pth = Path(p)
        if pth.is_file():
            if _matches_any_pattern(pth, patterns) and _matches_no_ignores(
                pth, ignore_patterns
            ):
                jobs = default_file_classifier(str(pth), plan, coord_target_panel)
                if work_dir:
                    for j in jobs:
                        j.context.add_metadata("work_dir", work_dir)
                # If this submission is fail-only BAMs, downstream errors are expected/noisy.
                if fail_only_bam_submission and _is_bam_path(pth):
                    for j in jobs:
                        j.context.add_metadata("fail_only_bam_submission", True)
                # Add reference genome to job metadata if available
                if coord_reference:
                    for j in jobs:
                        j.context.add_metadata("reference", coord_reference)
                # Add target panel to job metadata if available
                if coord_target_panel:
                    for j in jobs:
                        j.context.add_metadata("target_panel", coord_target_panel)
                seed_jobs += jobs
        elif pth.is_dir():
            # Stream directory walk in small batches; avoid building huge lists
            walker = pth.rglob("*") if recursive else pth.glob("*")
            batch: List[Job] = []
            for f in walker:
                if not f.is_file():
                    continue
                if not _matches_any_pattern(f, patterns) or not _matches_no_ignores(
                    f, ignore_patterns
                ):
                    continue
                jobs = default_file_classifier(str(f), plan, coord_target_panel)
                if work_dir:
                    for j in jobs:
                        j.context.add_metadata("work_dir", work_dir)
                if fail_only_bam_submission and _is_bam_path(f):
                    for j in jobs:
                        j.context.add_metadata("fail_only_bam_submission", True)
                # Add reference genome to job metadata if available
                if coord_reference:
                    for j in jobs:
                        j.context.add_metadata("reference", coord_reference)
                # Add target panel to job metadata if available
                if coord_target_panel:
                    for j in jobs:
                        j.context.add_metadata("target_panel", coord_target_panel)
                batch += jobs
                if len(batch) >= 256:
                    while True:
                        try:
                            can_accept = await coord.can_accept_jobs.remote(len(batch))
                        except Exception:
                            can_accept = True
                        if can_accept:
                            await coord.submit_jobs.remote(batch)
                            break
                        await asyncio.sleep(0.2)
                    batch = []
            if batch:
                while True:
                    try:
                        can_accept = await coord.can_accept_jobs.remote(len(batch))
                    except Exception:
                        can_accept = True
                    if can_accept:
                        await coord.submit_jobs.remote(batch)
                        break
                    await asyncio.sleep(0.2)
    if seed_jobs:
        # Submit in bounded batches to avoid flooding the coordinator/actors
        BATCH = 256
        for i in range(0, len(seed_jobs), BATCH):
            chunk = seed_jobs[i : i + BATCH]
            while True:
                try:
                    can_accept = await coord.can_accept_jobs.remote(len(chunk))
                except Exception:
                    can_accept = True
                if can_accept:
                    await coord.submit_jobs.remote(chunk)
                    break
                await asyncio.sleep(0.2)


async def tqdm_monitor(coord, continuous: bool = False) -> None:
    # Create queue-specific progress bars (similar to your original). If
    # continuous=True, do not exit when queues drain; keep running for
    # watchdog-driven new work until interrupted.
    order = [
        "preprocessing",
        "bed_conversion",
        "mgmt",
        "cnv",
        "target",
        "fusion",
        "classification",
        "slow",
    ]
    bars = {
        q: tqdm(desc=q.title(), unit="jobs", position=i, leave=True, ncols=120, 
                bar_format='{l_bar}{bar:20}{r_bar}')
        for i, q in enumerate(order)
    }
    overall = tqdm(
        desc="Overall Progress", unit="jobs", position=len(order), leave=True, ncols=120,
        bar_format='{l_bar}{bar:20}{r_bar}'
    )

    try:
        while True:
            s = await coord.stats.remote()
            # Build completed per queue (completed + failed)
            completed_by_type = s.get("completed_by_type", {})
            failed_by_type = s.get("failed_by_type", {})

            completed_queue_counts: Dict[str, int] = {q: 0 for q in order}
            for jt, c in completed_by_type.items():
                q = job_queue_of(jt)
                completed_queue_counts[q] += c
            for jt, c in failed_by_type.items():
                q = job_queue_of(jt)
                completed_queue_counts[q] += c

            # Active per queue
            active_by_queue = s.get("active_by_queue", {})
            active_count_by_queue = s.get("active_count_by_queue", {}) or {}
            active_count_by_job_type = s.get("active_count_by_job_type", {}) or {}
            waiting_count_by_queue = s.get("waiting_count_by_queue", {}) or {}
            active_count_by_queue = s.get("active_count_by_queue", {}) or {}
            active_count_by_job_type = s.get("active_count_by_job_type", {}) or {}
            waiting_count_by_queue = s.get("waiting_count_by_queue", {}) or {}

            # Update per-queue bars
            for q in order:
                completed_in_q = completed_queue_counts.get(q, 0)
                # Use full active counts for debugging; UI list is truncated to 2/jobs/queue.
                active_in_q = int(active_count_by_queue.get(q, 0) or 0)
                total_for_q = completed_in_q + active_in_q
                b = bars[q]
                b.total = max(b.total or 0, total_for_q)
                b.n = completed_in_q
                # Show up to 2 active jobs with per-job duration (no full path to avoid clutter)
                active_jobs_info = active_by_queue.get(q, [])
                if active_jobs_info:
                    parts = []
                    for j in active_jobs_info[:2]:
                        fn = os.path.basename(j.get("filepath", ""))
                        if len(fn) > 25:
                            fn = fn[:22] + "..."
                        dq = j.get("dispatch_wait_seconds")
                        if dq is not None and int(dq) > 0:
                            parts.append(f"{fn} ({j['duration']}s, +{int(dq)}q)")
                        else:
                            parts.append(f"{fn} ({j['duration']}s)")
                    b.set_postfix_str(" | ".join(parts))
                else:
                    b.set_postfix_str("")
                b.refresh()

            # Overall
            total_with_waiting = (
                s.get("total_enqueued", 0) - s.get("total_skipped", 0)
            ) + int(s.get("waiting_serialized", 0) or 0)
            overall.total = max(overall.total or 0, total_with_waiting) or 0
            overall.n = s["completed"] + s["failed"]
            fail_str = f"Fail:{s['failed']}"
            oom = s.get("oom_count", 0)
            if oom and s["failed"]:
                fail_str += f" (OOM:{oom})"
            top_active_types = sorted(
                active_count_by_job_type.items(), key=lambda x: x[1], reverse=True
            )[:3]
            top_waiting_queues = sorted(
                waiting_count_by_queue.items(), key=lambda x: x[1], reverse=True
            )[:3]
            active_types_str = ", ".join([f"{k}:{v}" for k, v in top_active_types])
            waiting_queues_str = ", ".join(
                [f"{k}:{v}" for k, v in top_waiting_queues]
            )
            overall.set_postfix_str(
                f"Act:{s['active_count']}[{active_types_str}] | "
                f"Wait:{s.get('waiting_global',0)}+{s.get('waiting_serialized',0)} "
                f"({waiting_queues_str}) | "
                f"Done:{s['completed']} | {fail_str} | "
                f"Skip:{s['total_skipped']} | Tot:{s['total_enqueued']}"
            )
            overall.refresh()

            # exit condition: only in non-continuous mode (batch). In continuous
            # mode, keep monitoring to allow watchdog to enqueue new files.
            if (
                (not continuous)
                and s["active_count"] == 0
                and overall.n >= overall.total
            ):
                break

            await asyncio.sleep(0.5)
    finally:
        for b in bars.values():
            b.close()
        overall.close()
        # Write failed jobs log when run ends so user can inspect (e.g. OOM)
        try:
            final = ray.get(coord.stats.remote())
            if final.get("failed", 0) > 0:
                path = ray.get(coord._write_failed_jobs_file.remote())
                if path:
                    oom = final.get("oom_count", 0)
                    print(f"\nFailed jobs log: {path} ({final['failed']} failures, {oom} likely OOM)")
        except Exception:
            pass


def _should_use_rich_progress() -> bool:
    if not _RICH_AVAILABLE:
        return False
    preference = os.environ.get("ROBIN_PROGRESS", "").strip().lower()
    if preference in {"tqdm", "plain", "off", "0", "false", "no"}:
        return False
    if preference in {"rich", "on", "1", "true", "yes"}:
        return True
    return True


async def rich_monitor(coord, continuous: bool = False) -> None:
    order = [
        "preprocessing",
        "bed_conversion",
        "mgmt",
        "cnv",
        "target",
        "fusion",
        "classification",
        "slow",
    ]
    progress = Progress(
        SpinnerColumn(),
        TextColumn("{task.description}"),
        BarColumn(bar_width=20),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        TextColumn("{task.fields[detail]}"),
        refresh_per_second=4,
    )
    task_ids = {
        q: progress.add_task(q.title(), total=None, detail="") for q in order
    }
    overall_id = progress.add_task("Overall Progress", total=None, detail="")

    if _RICH_CONSOLE:
        _RICH_CONSOLE.print(
            "[dim]Progress: time column = run elapsed; (Ns) = handler wall time once running "
            "(+Mq = seconds after dispatch waiting in Pool/Ray); before handler starts, wall time since dispatch[/dim]"
        )
    progress.start()
    try:
        while True:
            s = await coord.stats.remote()
            completed_by_type = s.get("completed_by_type", {})
            failed_by_type = s.get("failed_by_type", {})

            completed_queue_counts: Dict[str, int] = {q: 0 for q in order}
            for jt, c in completed_by_type.items():
                q = job_queue_of(jt)
                completed_queue_counts[q] += c
            for jt, c in failed_by_type.items():
                q = job_queue_of(jt)
                completed_queue_counts[q] += c

            active_by_queue = s.get("active_by_queue", {})
            active_count_by_queue = s.get("active_count_by_queue", {}) or {}
            active_count_by_job_type = s.get("active_count_by_job_type", {}) or {}
            waiting_count_by_queue = s.get("waiting_count_by_queue", {}) or {}

            for q in order:
                completed_in_q = completed_queue_counts.get(q, 0)
                # Use full active counts for debugging; UI list is truncated to 2/jobs/queue.
                active_in_q = int(active_count_by_queue.get(q, 0) or 0)
                total_for_q = completed_in_q + active_in_q
                active_jobs_info = active_by_queue.get(q, [])
                parts = []
                total_for_q_display = str(total_for_q) if total_for_q > 0 else "-"
                parts.append(f"{completed_in_q}/{total_for_q_display}")
                if active_jobs_info:
                    for j in active_jobs_info[:2]:
                        fn = os.path.basename(j.get("filepath", ""))
                        if len(fn) > 25:
                            fn = fn[:22] + "..."
                        dq = j.get("dispatch_wait_seconds")
                        if dq is not None and int(dq) > 0:
                            parts.append(f"{fn} ({j['duration']}s, +{int(dq)}q)")
                        else:
                            parts.append(f"{fn} ({j['duration']}s)")
                detail = " | ".join(parts)
                progress.update(
                    task_ids[q],
                    total=total_for_q if total_for_q > 0 else None,
                    completed=completed_in_q,
                    detail=detail,
                )

            total_with_waiting = (
                s.get("total_enqueued", 0) - s.get("total_skipped", 0)
            ) + int(s.get("waiting_serialized", 0) or 0)
            fail_detail = f"Fail:{s['failed']}"
            oom = s.get("oom_count", 0)
            if oom and s.get("failed"):
                fail_detail += f" (OOM:{oom})"
            progress.update(
                overall_id,
                total=total_with_waiting if total_with_waiting > 0 else None,
                completed=s["completed"] + s["failed"],
                detail=(
                    f"Done:{s['completed'] + s['failed']}/{total_with_waiting} | "
                    f"Act:{s['active_count']}[{', '.join([f'{k}:{v}' for k, v in sorted(active_count_by_job_type.items(), key=lambda x: x[1], reverse=True)[:3]])}] | "
                    f"Wait:{s.get('waiting_global',0)}+{s.get('waiting_serialized',0)} "
                    f"({', '.join([f'{k}:{v}' for k, v in sorted(waiting_count_by_queue.items(), key=lambda x: x[1], reverse=True)[:3]])}) | "
                    f"{fail_detail} | Skip:{s['total_skipped']}"
                ),
            )

            if (
                (not continuous)
                and s["active_count"] == 0
                and s["completed"] + s["failed"] >= total_with_waiting
            ):
                break

            await asyncio.sleep(0.5)
    finally:
        progress.stop()
        # Write failed jobs log when run ends so user can inspect (e.g. OOM)
        try:
            final = ray.get(coord.stats.remote())
            if final.get("failed", 0) > 0:
                path = ray.get(coord._write_failed_jobs_file.remote())
                if path:
                    oom = final.get("oom_count", 0)
                    print(f"\nFailed jobs log: {path} ({final['failed']} failures, {oom} likely OOM)")
        except Exception:
            pass


class RayFileWatcher(FileSystemEventHandler):
    def __init__(
        self,
        coord,
        plan: List[str],
        target_panel: str,
        patterns: Optional[List[str]] = None,
        ignore_patterns: Optional[List[str]] = None,
        recursive: bool = True,
        work_dir: Optional[str] = None,
    ):
        self.coord = coord
        self.plan = plan
        self.target_panel = target_panel
        self.patterns = patterns or ["*"]
        self.ignore_patterns = ignore_patterns or []
        self.recursive = recursive
        self.work_dir = work_dir
        self.processed: set[str] = set()
        # Rate limiting / batching
        self._pending_jobs: List[Job] = []
        self._last_flush: float = time.time()
        self._flush_interval_s: float = 0.5
        self._batch_size: int = 128
        # Track whether this watch session has seen any "pass" BAMs.
        # If a pass BAM is present, we want fail-BAM errors to remain visible.
        self._watch_seen_bam: bool = False
        self._watch_seen_pass_bam: bool = False

    def _should_process(self, fp: str) -> bool:
        p = Path(fp)
        if any(p.match(ip) for ip in self.ignore_patterns):
            return False
        if self.patterns:
            return any(p.match(pt) for pt in self.patterns)
        return True

    def _flush_if_needed(self, force: bool = False) -> None:
        try:
            if not self._pending_jobs:
                return
            now = time.time()
            if (
                force
                or len(self._pending_jobs) >= self._batch_size
                or (now - self._last_flush) >= self._flush_interval_s
            ):
                batch = self._pending_jobs[: self._batch_size]
                del self._pending_jobs[: self._batch_size]
                try:
                    ray.get(self.coord.submit_jobs.remote(batch))
                except Exception:
                    pass
                self._last_flush = now
        except Exception:
            pass

    def _handle(self, fp: str):
        if fp in self.processed:
            return
        if not self._should_process(fp):
            return
        self.processed.add(fp)
        # Maintain watch-level pass/fail BAM state
        try:
            if _is_bam_path(Path(fp)):
                self._watch_seen_bam = True
                k = _bam_filename_kind(fp)
                if k == "pass":
                    self._watch_seen_pass_bam = True
        except Exception:
            pass
        jobs = default_file_classifier(fp, self.plan, self.target_panel)
        if self.work_dir:
            for j in jobs:
                j.context.add_metadata("work_dir", self.work_dir)
        # Tag fail-only BAM submissions for this watch session.
        # If we've seen a pass BAM, do not suppress errors.
        try:
            if _is_bam_path(Path(fp)) and self._watch_seen_bam and (not self._watch_seen_pass_bam):
                if _bam_filename_kind(fp) == "fail":
                    for j in jobs:
                        j.context.add_metadata("fail_only_bam_submission", True)
        except Exception:
            pass
        # enqueue and flush under rate limiter
        self._pending_jobs.extend(jobs)
        self._flush_if_needed()

    def on_created(self, event):
        if not event.is_directory:
            self._handle(event.src_path)

    def on_modified(self, event):
        if not event.is_directory:
            self._handle(event.src_path)

    def on_moved(self, event):
        """Handle file move events (common with rsync default behavior)."""
        if not event.is_directory:
            self._handle(event.dest_path)  # Use dest_path for the final location

    # Ensure remaining jobs are flushed on shutdown
    def flush_remaining(self):
        while self._pending_jobs:
            self._flush_if_needed(force=True)


# Global coordinator reference for external access
_GLOBAL_COORDINATOR = None

# Global observer and watcher for adding paths at runtime (GUI "Add Folder to Watch")
_GLOBAL_OBSERVER = None
_GLOBAL_WATCHER = None
_GLOBAL_WATCH_CONTEXT: Optional[Dict[str, Any]] = None
# Maps normalized path string -> ObservedWatch for add/remove
_GLOBAL_WATCHED_PATHS: Dict[str, Any] = {}


def get_watched_paths() -> List[str]:
    """Return list of currently watched directory paths (normalized)."""
    global _GLOBAL_WATCHED_PATHS
    return sorted(_GLOBAL_WATCHED_PATHS.keys())


def remove_watch_path(path_to_remove: str) -> Tuple[bool, str]:
    """
    Remove a directory from the workflow's watch list.

    Returns:
        Tuple of (success: bool, message: str)
    """
    global _GLOBAL_OBSERVER, _GLOBAL_WATCHED_PATHS

    if _GLOBAL_OBSERVER is None:
        return False, "No workflow observer found. Is a Ray workflow with watch enabled running?"

    pth = Path(path_to_remove).resolve()
    path_key = str(pth)

    if path_key not in _GLOBAL_WATCHED_PATHS:
        return False, f"Path is not being watched: {path_to_remove}"

    watch = _GLOBAL_WATCHED_PATHS[path_key]
    try:
        _GLOBAL_OBSERVER.unschedule(watch)
        del _GLOBAL_WATCHED_PATHS[path_key]
        return True, f"Removed watch path: {path_to_remove}"
    except Exception as e:
        return False, f"Failed to remove watch: {e}"


def _path_overlaps_work_dir(candidate: Path, work_dir: Optional[str]) -> bool:
    """
    Return True if candidate path overlaps with work directory (should reject).

    Rejects when:
    - candidate equals work_dir
    - candidate is inside work_dir (watching our own output)
    - work_dir is inside candidate (would watch work_dir and its outputs)
    """
    if not work_dir:
        return False
    try:
        wd = Path(work_dir).resolve()
        cand = candidate.resolve()
        if cand == wd:
            return True
        # Reject if candidate is inside work_dir (watching our own output)
        cand.relative_to(wd)
        return True
    except ValueError:
        pass
    try:
        wd = Path(work_dir).resolve()
        cand = candidate.resolve()
        # Reject if work_dir is inside candidate (would watch work_dir)
        wd.relative_to(cand)
        return True
    except ValueError:
        pass
    return False


def _scan_watch_folder_for_sample_dirs(
    folder: Path,
    work_dir: Optional[str],
    patterns: Optional[List[str]],
    ignore_patterns: Optional[List[str]],
    recursive: bool,
    max_bams_to_probe: int = 5000,
    max_bams_per_dir: int = 8,
) -> Tuple[List[Path], List[str], Optional[str]]:
    """
    Scan a candidate watch folder and return:
    - A list of directory paths that are safe to watch (contain BAMs for samples
      that do NOT appear previously analysed, and contain no BAMs for analysed samples)
    - A sorted list of analysed sample IDs detected
    - Optional error message

    Returns:
        (dirs_to_watch, analysed_sample_ids, error_message)
    """
    if not work_dir:
        return [], [], None

    try:
        wd = Path(work_dir).resolve()
    except Exception:
        return [], [], None

    if not wd.exists() or not wd.is_dir():
        return [], [], None

    def _sample_dir_has_analysis_artifacts_local(sample_dir: Path) -> bool:
        """Heuristic: True if sample_dir contains real analysis outputs."""
        try:
            if not sample_dir.exists() or not sample_dir.is_dir():
                return False

            # Strong indicators (target outputs / SNP inputs)
            if (sample_dir / "target.bam").exists() or (sample_dir / "target.bam.bai").exists():
                return True
            if (sample_dir / "targets_exceeding_threshold.bed").exists():
                return True

            # Fusion outputs
            if (sample_dir / "target_fusion.csv").exists() or (sample_dir / "genome_wide_fusion.csv").exists():
                return True

            # Classifier outputs
            if (sample_dir / "sturgeon_scores.csv").exists():
                return True
            if (sample_dir / "NanoDX_scores.csv").exists() or (sample_dir / "PanNanoDX_scores.csv").exists():
                return True
            if (sample_dir / "random_forest_scores.csv").exists():
                return True

            # Common analysis directories
            if (sample_dir / "bed_files").is_dir():
                return True
            if (sample_dir / "igv").is_dir():
                return True
            if (sample_dir / "clair3").is_dir():
                return True

            # CNV outputs
            try:
                for child in sample_dir.glob("*_copy_numbers.pkl"):
                    if child.is_file():
                        return True
            except Exception:
                pass

            return False
        except Exception:
            return False

    # Import lazily so workflow can run even if analysis deps aren't present.
    try:
        from robin.analysis.bam_preprocessor import _extract_sample_id_from_bam  # type: ignore
    except Exception:
        return [], [], None

    analysed_sample_ids: Set[str] = set()
    unanalysed_sample_ids: Set[str] = set()
    # Map directory containing BAMs -> set(sample_ids present in that directory)
    dir_samples: Dict[Path, Set[str]] = {}
    # To keep folder-add responsive, do not probe every BAM in a directory.
    # We sample up to max_bams_per_dir BAMs per parent directory. If a directory
    # contains multiple samples mixed together, this will usually be detected quickly
    # and we will avoid auto-watching that directory.
    dir_probe_counts: Dict[Path, int] = {}
    bams_probed = 0
    try:
        walker = folder.rglob("*") if recursive else folder.glob("*")
        for f in walker:
            if not f.is_file():
                continue
            if not _matches_any_pattern(f, patterns) or not _matches_no_ignores(
                f, ignore_patterns
            ):
                continue

            try:
                parent_dir = f.parent.resolve()
            except Exception:
                parent_dir = f.parent

            # Limit per-directory probing to keep scan fast.
            seen = int(dir_probe_counts.get(parent_dir, 0))
            if seen >= max_bams_per_dir:
                continue

            bams_probed += 1
            if bams_probed > max_bams_to_probe:
                break

            try:
                sid = _extract_sample_id_from_bam(str(f))
            except Exception:
                sid = "unknown"

            if not sid or sid == "unknown":
                dir_probe_counts[parent_dir] = seen + 1
                continue

            # Track which sample IDs appear in which directories.
            dir_samples.setdefault(parent_dir, set()).add(sid)
            dir_probe_counts[parent_dir] = seen + 1

            # Classify sample_id by presence of analysis artifacts in work_dir
            sample_dir = wd / sid
            if sample_dir.is_dir() and _sample_dir_has_analysis_artifacts_local(sample_dir):
                analysed_sample_ids.add(sid)
            else:
                unanalysed_sample_ids.add(sid)
    except Exception as e:
        return [], [], f"{e}"

    analysed_sorted = sorted(analysed_sample_ids)

    # Choose directories to watch:
    # Only watch directories that contain BAMs exclusively for unanalysed samples.
    # This prevents reprocessing if a directory contains mixed samples.
    dirs_to_watch: List[Path] = []
    for d, sids_in_dir in dir_samples.items():
        if not sids_in_dir:
            continue
        if any((sid in analysed_sample_ids) for sid in sids_in_dir):
            continue
        # Must contain at least one unanalysed sample
        if any((sid in unanalysed_sample_ids) for sid in sids_in_dir):
            dirs_to_watch.append(d)

    # Stable ordering for messages/scheduling.
    dirs_to_watch = sorted(set(dirs_to_watch), key=lambda p: str(p))
    return dirs_to_watch, analysed_sorted, None


def add_watch_path(new_path: str) -> Tuple[bool, str]:
    """
    Add a directory to the running workflow's watch list.

    Schedules the path for file watching and submits existing BAM files for processing.
    Uses the same plan, work_dir, panel, and patterns as the original workflow.

    Returns:
        Tuple of (success: bool, message: str)
    """
    global _GLOBAL_COORDINATOR, _GLOBAL_OBSERVER, _GLOBAL_WATCHER, _GLOBAL_WATCH_CONTEXT, _GLOBAL_WATCHED_PATHS

    coord = get_coordinator_sync()
    if coord is None:
        return False, "No workflow coordinator found. Is a Ray workflow running?"

    if _GLOBAL_WATCH_CONTEXT is None:
        return False, "Workflow watch context not available (watch may be disabled)."

    pth = Path(new_path).resolve()
    if not pth.exists():
        return False, f"Path does not exist: {new_path}"
    if not pth.is_dir():
        return False, f"Path is not a directory: {new_path}"

    ctx = _GLOBAL_WATCH_CONTEXT
    plan = ctx.get("plan") or []
    work_dir = ctx.get("work_dir")

    # Reject paths that overlap with work directory (would watch already-analysed output)
    if _path_overlaps_work_dir(pth, work_dir):
        return False, (
            f"Cannot watch '{new_path}': it overlaps with the work directory. "
            "Watched folders must not contain or be inside the analysis output directory."
        )
    patterns = ctx.get("patterns") or ["*.bam"]
    ignore_patterns = ctx.get("ignore_patterns") or []
    recursive = ctx.get("recursive", True)

    # Option B behavior:
    # When a user adds a folder that contains a mix of previously-analysed and new samples,
    # automatically add only those subfolders that contain exclusively un-analysed samples.
    dirs_to_watch, analysed_sample_ids, scan_err = _scan_watch_folder_for_sample_dirs(
        folder=pth,
        work_dir=work_dir,
        patterns=patterns,
        ignore_patterns=ignore_patterns,
        recursive=recursive,
    )
    if scan_err:
        return False, (
            f"Cannot watch '{new_path}': error while scanning for existing analysed samples: {scan_err}"
        )
    if analysed_sample_ids and not dirs_to_watch:
        shown = analysed_sample_ids[:10]
        more = len(analysed_sample_ids) - len(shown)
        shown_str = ", ".join(shown)
        more_str = f" (+{more} more)" if more > 0 else ""
        return False, (
            f"Cannot watch '{new_path}': it contains only previously-analysed sample(s) (no new samples found to watch). "
            f"Samples: {shown_str}{more_str}. "
            "Move or delete the corresponding sample output folder(s) in the work directory if you want to reanalyse, "
            "or choose a different watch folder."
        )

    # If there are no BAMs (or no sample IDs extracted), fall back to watching the path as-is.
    if not analysed_sample_ids and not dirs_to_watch:
        dirs_to_watch = [pth]

    # Schedule for watching if observer is active
    newly_added: List[Path] = []
    already_watched: List[Path] = []
    if _GLOBAL_OBSERVER is not None and _GLOBAL_WATCHER is not None:
        for d in dirs_to_watch:
            path_key = str(Path(d).resolve())
            if path_key in _GLOBAL_WATCHED_PATHS:
                already_watched.append(d)
                continue
            try:
                watch = _GLOBAL_OBSERVER.schedule(
                    _GLOBAL_WATCHER, str(d), recursive=recursive
                )
                _GLOBAL_WATCHED_PATHS[path_key] = watch
                newly_added.append(d)
            except Exception as e:
                return False, f"Failed to schedule watch for '{d}': {e}"

    # Submit existing paths for processing (run in a thread to avoid event loop conflict
    # with NiceGUI or other async frameworks already running an event loop)
    _submit_error: List[Optional[Exception]] = [None]

    def _run_submit() -> None:
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        try:
            loop.run_until_complete(
                submit_existing_paths(
                    coord,
                    [str(d) for d in dirs_to_watch],
                    plan,
                    patterns=patterns,
                    ignore_patterns=ignore_patterns,
                    recursive=recursive,
                    work_dir=work_dir,
                )
            )
        except Exception as e:
            _submit_error[0] = e
        finally:
            loop.close()

    thread = threading.Thread(target=_run_submit)
    thread.start()
    thread.join()

    if _submit_error[0] is not None:
        return False, f"Failed to submit existing files: {_submit_error[0]}"

    # Build a user-facing message
    msg_parts: List[str] = []
    if newly_added:
        if len(newly_added) <= 5:
            msg_parts.append("Added watch path(s): " + "; ".join(str(p) for p in newly_added))
        else:
            msg_parts.append(f"Added {len(newly_added)} watch path(s) under: {new_path}")
    if already_watched:
        msg_parts.append(f"{len(already_watched)} path(s) were already watched")
    if analysed_sample_ids:
        shown = analysed_sample_ids[:10]
        more = len(analysed_sample_ids) - len(shown)
        shown_str = ", ".join(shown)
        more_str = f" (+{more} more)" if more > 0 else ""
        msg_parts.append(
            f"Skipped previously-analysed sample(s): {shown_str}{more_str} (not added to watch)"
        )
    if not msg_parts:
        msg_parts.append(f"Added watch path: {new_path}")

    return True, " | ".join(msg_parts)


async def _get_coordinator():
    """Get the global coordinator reference for external job submission."""
    global _GLOBAL_COORDINATOR
    if _GLOBAL_COORDINATOR is None:
        # Try to get the coordinator by name if it exists
        try:
            _GLOBAL_COORDINATOR = ray.get_actor("robin_coordinator")
        except Exception:
            pass
    return _GLOBAL_COORDINATOR


def get_coordinator_sync():
    """Get the global coordinator reference synchronously."""
    global _GLOBAL_COORDINATOR
    if _GLOBAL_COORDINATOR is None:
        # Try to get the coordinator by name if it exists
        try:
            _GLOBAL_COORDINATOR = ray.get_actor("robin_coordinator")
        except Exception:
            pass
    return _GLOBAL_COORDINATOR


def set_coordinator(coord):
    """Set the global coordinator reference from outside."""
    global _GLOBAL_COORDINATOR
    _GLOBAL_COORDINATOR = coord


async def run(
    plan: List[str],
    paths: List[str],
    target_panel: str,
    analysis_workers: int = 1,
    process_existing: bool = True,
    monitor: bool = True,
    watch: bool = False,
    patterns: Optional[List[str]] = None,
    ignore_patterns: Optional[List[str]] = None,
    recursive: bool = True,
    work_dir: Optional[str] = None,
    log_level: str = "INFO",
    preset: Optional[str] = None,
    workflow_runner: Any = None,
    reference: Optional[str] = None,
    gui_host: str = "0.0.0.0",
    gui_port: int = 8081,
    center: str = None,
    enable_batching: bool = True,
):
    global GLOBAL_LOG_LEVEL, _GLOBAL_OBSERVER, _GLOBAL_WATCHER, _GLOBAL_WATCH_CONTEXT, _GLOBAL_WATCHED_PATHS
    GLOBAL_LOG_LEVEL = (log_level or "INFO").upper()

    # Configure Ray logging to reduce verbose output
    import logging

    ray_logger = logging.getLogger("ray")
    ray_logger.setLevel(logging.WARNING)

    # Also reduce Raylet logging
    raylet_logger = logging.getLogger("raylet")
    raylet_logger.setLevel(logging.WARNING)

    # Reduce other Ray-related logging
    logging.getLogger("ray.worker").setLevel(logging.WARNING)
    logging.getLogger("ray.remote").setLevel(logging.WARNING)
    logging.getLogger("ray.actor").setLevel(logging.WARNING)
    logging.getLogger("ray.util").setLevel(logging.WARNING)

    # Set Ray environment variables to reduce verbose output
    import os

    os.environ["RAY_DISABLE_IMPORT_WARNING"] = "1"
    os.environ["RAY_DISABLE_DEPRECATION_WARNING"] = "1"

    # Reference genome status (minimal logging)
    if reference:
        pass  # Reference genome provided
    else:
        pass  # No reference genome provided

    # Ensure any previous coordinator is terminated to avoid stale-code actors
    try:
        old = ray.get_actor("robin_coordinator")
        try:
            ray.kill(old)
        except Exception:
            pass
    except Exception:
        pass
    # Start a fresh coordinator (not detached)
    job_timeout = _job_timeout_seconds()
    if job_timeout > 0:
        logging.info(f"Per-job timeout enabled: jobs running longer than {job_timeout}s will be cancelled (ROBIN_JOB_TIMEOUT_SECONDS)")
    try:
        coord = Coordinator.options(name="robin_coordinator", num_cpus=0).remote(
            target_panel=target_panel, analysis_workers=analysis_workers, preset=preset, reference=reference, enable_batching=enable_batching, job_timeout_seconds=job_timeout
        )
    except Exception:
        coord = Coordinator.options(name="robin_coordinator").remote(
            target_panel=target_panel, analysis_workers=analysis_workers, preset=preset, reference=reference, enable_batching=enable_batching, job_timeout_seconds=job_timeout
        )

    # Set global coordinator reference for external access
    global _GLOBAL_COORDINATOR
    _GLOBAL_COORDINATOR = coord

    # run async setup without blocking the event loop
    try:
        await coord.setup.remote()
    except Exception:
        pass
    
    # Set work_dir on coordinator so it's available during shutdown
    if work_dir:
        try:
            coord.set_work_dir.remote(work_dir)
        except Exception:
            pass

    # Launch GUI early - immediately after coordinator setup but before job submission
    # This ensures the GUI is available as soon as possible, before any jobs start
    gui_launcher = None
    gui_publish_task = None
    try:
        if _GUIUpdateType is None:
            print("GUI not launched: GUI modules unavailable on this environment.")
        elif not work_dir:
            print("GUI not launched: --work-dir not provided.")
        else:
            _print_styled(
                "Launching GUI early to ensure immediate availability...", level="info"
            )
            launcher = _gui_launch(
                host=gui_host,
                port=gui_port,
                workflow_runner=workflow_runner,
                workflow_steps=plan,
                monitored_directory=work_dir,
                center=center,
            )
            gui_launcher = launcher
            try:
                url = (
                    launcher.get_gui_url()
                    if hasattr(launcher, "get_gui_url")
                    else "http://localhost:8081"
                )
                _print_styled(f"GUI launched successfully on {url}", level="success")
            except Exception:
                _print_styled("GUI launched successfully", level="success")

            # Start GUI update publishing task immediately
            async def _publish_gui():
                while True:
                    try:
                        s = await coord.stats.remote()
                        # Basic workflow status & progress
                        total = int(
                            (s.get("total_enqueued", 0) - s.get("total_skipped", 0))
                            or 0
                        )
                        completed = int(s.get("completed", 0) or 0)
                        failed = int(s.get("failed", 0) or 0)
                        progress = (completed + failed) / total if total > 0 else 0.0
                        _gui_send_update(
                            _GUIUpdateType.WORKFLOW_STATUS,
                            {"is_running": True, "start_time": time.time()},
                            priority=1,
                        )
                        _gui_send_update(
                            _GUIUpdateType.PROGRESS_UPDATE,
                            {
                                "progress": progress,
                                "completed": completed,
                                "failed": failed,
                                "total": total,
                            },
                            priority=1,
                        )

                        # Queue status (running/total by category). If the coordinator
                        # doesn't provide category totals/running, derive them from the
                        # available fields so the GUI can still display correct values.
                        ru = s.get("running_by_category", {}) or {}
                        tu = s.get("totals_by_category", {}) or {}

                        # Derive running by category when missing using active_by_queue
                        if not ru or sum(int(v or 0) for v in ru.values()) == 0:
                            ru = {
                                "preprocessing": 0,
                                "mgmt": 0,
                                "cnv": 0,
                                "target": 0,
                                "fusion": 0,
                                "classification": 0,
                                "other": 0,
                            }
                            abq = s.get("active_by_queue", {}) or {}
                            for qname, jobs in abq.items():
                                n = len(jobs) if isinstance(jobs, list) else 0
                                if qname in {
                                    "preprocessing",
                                    "mgmt",
                                    "cnv",
                                    "target",
                                    "fusion",
                                    "classification",
                                }:
                                    ru[qname] += n
                                elif qname in {"bed_conversion", "slow", "other"}:
                                    ru["other"] += n

                        # Derive totals by category when missing using completed/failed + active
                        if not tu or sum(int(v or 0) for v in tu.values()) == 0:
                            tu = {
                                "preprocessing": 0,
                                "mgmt": 0,
                                "cnv": 0,
                                "target": 0,
                                "fusion": 0,
                                "classification": 0,
                                "other": 0,
                            }
                            cbt = s.get("completed_by_type", {}) or {}
                            fbt = s.get("failed_by_type", {}) or {}

                            def _cat_of_local(jt: str) -> str:
                                if jt == "preprocessing":
                                    return "preprocessing"
                                if jt in {"mgmt", "cnv", "target", "fusion"}:
                                    return jt
                                if jt in CLASSIFICATION_TYPES:
                                    return "classification"
                                return "other"

                            for jt, c in cbt.items():
                                tu[_cat_of_local(jt)] += int(c or 0)
                            for jt, c in fbt.items():
                                tu[_cat_of_local(jt)] += int(c or 0)
                            # Add currently active jobs
                            abq2 = s.get("active_by_queue", {}) or {}
                            for qname, jobs in abq2.items():
                                n = len(jobs) if isinstance(jobs, list) else 0
                                if qname in {
                                    "preprocessing",
                                    "mgmt",
                                    "cnv",
                                    "target",
                                    "fusion",
                                    "classification",
                                }:
                                    tu[qname] += n
                                elif qname in {"bed_conversion", "slow", "other"}:
                                    tu["other"] += n

                        queue_payload = {
                            "preprocessing": {
                                "running": int(ru.get("preprocessing", 0) or 0),
                                "total": int(tu.get("preprocessing", 0) or 0),
                            },
                            "mgmt": {
                                "running": int(ru.get("mgmt", 0) or 0),
                                "total": int(tu.get("mgmt", 0) or 0),
                            },
                            "cnv": {
                                "running": int(ru.get("cnv", 0) or 0),
                                "total": int(tu.get("cnv", 0) or 0),
                            },
                            "target": {
                                "running": int(ru.get("target", 0) or 0),
                                "total": int(tu.get("target", 0) or 0),
                            },
                            "fusion": {
                                "running": int(ru.get("fusion", 0) or 0),
                                "total": int(tu.get("fusion", 0) or 0),
                            },
                            "classification": {
                                "running": int(ru.get("classification", 0) or 0),
                                "total": int(tu.get("classification", 0) or 0),
                            },
                            "other": {
                                "running": int(ru.get("other", 0) or 0),
                                "total": int(tu.get("other", 0) or 0),
                            },
                        }
                        _gui_send_update(
                            _GUIUpdateType.QUEUE_UPDATE, queue_payload, priority=2
                        )
                        # Also emit a human-readable log line so users see queue summaries in Live Logs
                        try:
                            pre = queue_payload["preprocessing"]
                            an_run = sum(
                                int(queue_payload[q]["running"])
                                for q in ["mgmt", "cnv", "target", "fusion"]
                            )
                            an_tot = sum(
                                int(queue_payload[q]["total"])
                                for q in ["mgmt", "cnv", "target", "fusion"]
                            )
                            cl = queue_payload["classification"]
                            ot = queue_payload["other"]
                            summary = (
                                f"Queues | Pre:{pre['running']}/{pre['total']} "
                                f"| An:{an_run}/{an_tot} "
                                f"| Cl:{cl['running']}/{cl['total']} "
                                f"| Ot:{ot['running']}/{ot['total']}"
                            )
                            _gui_send_update(
                                _GUIUpdateType.LOG_MESSAGE,
                                {"message": summary, "level": "INFO"},
                                priority=0,
                            )
                        except Exception:
                            pass

                        # Active jobs table
                        rows = []
                        for qname, jobs in (s.get("active_by_queue", {}) or {}).items():
                            for j in jobs:
                                rj = {
                                    "job_id": str(j.get("job_id", "")),
                                    "job_type": j.get("job_type", ""),
                                    "filepath": j.get("filepath", ""),
                                    "worker_name": qname,
                                    "duration": int(j.get("duration", 0) or 0),
                                    "progress": 0.0,
                                }
                                dw = j.get("dispatch_wait_seconds")
                                if dw is not None:
                                    rj["dispatch_wait_seconds"] = int(dw)
                                rows.append(rj)
                        _gui_send_update(
                            _GUIUpdateType.JOB_UPDATE, {"active_jobs": rows}, priority=1
                        )

                        # Samples overview
                        samples = s.get("samples", []) or []
                        _gui_send_update(
                            _GUIUpdateType.SAMPLES_UPDATE,
                            {"samples": samples},
                            priority=1,
                        )

                        await asyncio.sleep(15.0)  # Poll every 15 seconds to reduce update frequency
                    except asyncio.CancelledError:
                        break
                    except Exception:
                        await asyncio.sleep(15.0)  # Poll every 15 seconds to reduce update frequency

            gui_publish_task = asyncio.create_task(_publish_gui())
            _print_styled(
                "GUI monitoring started - workflow status will be updated in real-time",
                level="info",
            )
    except Exception as e:
        print(f"Warning: GUI failed to launch: {e}")

    # Create monitoring tasks BEFORE submitting paths so progress appears immediately
    tasks = []
    if monitor:
        monitor_fn = rich_monitor if _should_use_rich_progress() else tqdm_monitor
        tasks.append(asyncio.create_task(monitor_fn(coord, continuous=watch)))
    
    # Add GUI publish task if GUI was launched
    if gui_publish_task is not None:
        tasks.append(gui_publish_task)

    # Now submit existing paths (this will run while monitor displays progress)
    if process_existing and paths:
        await submit_existing_paths(
            coord,
            paths,
            plan,
            patterns=patterns,
            ignore_patterns=ignore_patterns,
            recursive=recursive,
            work_dir=work_dir,
        )

    observer = None
    watcher = None
    if watch and paths:
        observer = Observer()
        watcher = RayFileWatcher(
            coord,
            plan,
            target_panel,
            patterns=patterns,
            ignore_patterns=ignore_patterns,
            recursive=recursive,
            work_dir=work_dir,
        )
        for p in paths:
            if Path(p).is_dir():
                watch = observer.schedule(watcher, p, recursive=recursive)
                _GLOBAL_WATCHED_PATHS[str(Path(p).resolve())] = watch
        observer.start()

        # Store globally for GUI "Add Folder to Watch"
        _GLOBAL_OBSERVER = observer
        _GLOBAL_WATCHER = watcher
        _GLOBAL_WATCH_CONTEXT = {
            "plan": plan,
            "work_dir": work_dir,
            "patterns": patterns,
            "ignore_patterns": ignore_patterns,
            "recursive": recursive,
        }

    try:
        if tasks:
            await asyncio.gather(*tasks)
    except KeyboardInterrupt:
        print("\n[SHUTDOWN] Interrupted by user (Ctrl-C)")
        print("[SHUTDOWN] Initiating graceful shutdown...")
        # Cancel all running tasks gracefully
        if tasks:
            running_tasks = [t for t in tasks if not t.done()]
            if running_tasks:
                print(f"[SHUTDOWN] Cancelling {len(running_tasks)} running task(s)...")
                for task in running_tasks:
                    task.cancel()
                    try:
                        await task
                    except asyncio.CancelledError:
                        pass
                print("[SHUTDOWN] Tasks cancelled")
        # Signal coordinator to stop accepting new jobs and finalize
        try:
            print("[SHUTDOWN] Shutting down coordinator...")
            await coord.shutdown.remote()
        except Exception as e:
            print(f"[SHUTDOWN] Warning: Error shutting down coordinator: {e}")
        # Kill all Ray actors gracefully
        try:
            ray.kill(coord)
        except Exception:
            pass
        print("[SHUTDOWN] Workflow stopped gracefully")
    finally:
        # Always finalize target BAMs on exit (including normal completion)
        try:
            await coord.shutdown.remote()
        except Exception:
            pass
        if watcher is not None:
            try:
                print("[SHUTDOWN] Flushing file watcher...")
                watcher.flush_remaining()
                print("[SHUTDOWN] File watcher flushed")
            except Exception as e:
                print(f"[SHUTDOWN] Warning: Error flushing file watcher: {e}")
        if observer is not None:
            print("[SHUTDOWN] Stopping file observer...")
            observer.stop()
            observer.join()
            print("[SHUTDOWN] File observer stopped")

        # Clear global watch state
        _GLOBAL_OBSERVER = None
        _GLOBAL_WATCHER = None
        _GLOBAL_WATCH_CONTEXT = None
        _GLOBAL_WATCHED_PATHS.clear()

        # Final cleanup: ensure Ray is properly shut down
        try:
            import ray
            if ray.is_initialized():
                print("[SHUTDOWN] Shutting down Ray...")
                ray.shutdown()
                print("[SHUTDOWN] Ray shutdown complete")
        except Exception as e:
            print(f"[SHUTDOWN] Warning: Error shutting down Ray: {e}")
        print("[SHUTDOWN] All cleanup operations completed")


def parse_args():
    p = argparse.ArgumentParser(description="Ray Core robin workflow")
    p.add_argument(
        "--plan",
        action="append",
        default=[],
        help="Workflow plan entries like 'queue:job_type' in order",
    )
    p.add_argument(
        "--paths",
        nargs="*",
        default=[],
        help="Files/dirs to process immediately (optional)",
    )
    p.add_argument(
        "--analysis-workers",
        type=int,
        default=1,
        help="Parallel workers per analysis queue",
    )
    p.add_argument("--no-monitor", action="store_true", help="Disable tqdm monitor")
    p.add_argument(
        "--no-process-existing",
        action="store_true",
        help="Do not process existing files/dirs",
    )
    p.add_argument(
        "--watch",
        action="store_true",
        help="Watch directories for new files (Watchdog)",
    )
    p.add_argument(
        "--patterns",
        action="append",
        default=None,
        help="Glob patterns to include (repeatable)",
    )
    p.add_argument(
        "--ignore",
        action="append",
        default=None,
        help="Glob patterns to ignore (repeatable)",
    )
    p.add_argument(
        "--no-recursive",
        action="store_true",
        help="Disable recursive directory watching",
    )
    p.add_argument(
        "--ray-address", default=None, help="Ray cluster address (None = local)"
    )
    p.add_argument(
        "--no-clinvar-download",
        action="store_true",
        help="Do not auto-download ClinVar VCF files if missing",
    )
    p.add_argument(
        "--update-clinvar-if-newer",
        action="store_true",
        help="If available, update ClinVar to the newest NCBI version (best-effort Last-Modified check)",
    )
    p.add_argument(
        "--preset",
        default=None,
        choices=["p2i", "standard", "high"],
        help="Execution preset controlling actor grouping and concurrency",
    )
    p.add_argument(
        "--no-ray-dashboard",
        action="store_true",
        help="Disable Ray dashboard (default: dashboard enabled)",
    )
    p.add_argument(
        "--center",
        required=True,
        help="Center ID running the analysis (e.g., 'Oxford', 'Cambridge', 'London')",
    )
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Ensure ClinVar is present for variant annotation/reporting.
    if not args.no_clinvar_download:
        from robin.utils.clinvar_manager import (
            ensure_clinvar_files,
            update_clinvar_if_newer,
        )

        ensure_clinvar_files(download_if_missing=True)
        if args.update_clinvar_if_newer:
            update_clinvar_if_newer(download_if_missing=True)

    # Raise dashboard/State API job list limit (Ray default 10k) if not set
    if os.environ.get("RAY_MAX_LIMIT_FROM_DATA_SOURCE") is None:
        os.environ["RAY_MAX_LIMIT_FROM_DATA_SOURCE"] = "100000"
    try:
        # Expose dashboard on all interfaces when supported
        init_kwargs = {"address": args.ray_address}

        if not args.no_ray_dashboard:
            init_kwargs.update(
                {
                    "include_dashboard": True,
                    "dashboard_host": os.environ.get("RAY_DASHBOARD_HOST", "0.0.0.0"),
                }
            )

        ray.init(**init_kwargs)
    except TypeError:
        # Older Ray versions may not support dashboard args
        ray.init(address=args.ray_address)

    if not args.plan:
        # default example
        args.plan = [
            "preprocessing:preprocessing",
            "bed_conversion:bed_conversion",
            "mgmt:mgmt",
            "sturgeon:sturgeon",
        ]

    try:
        asyncio.run(
            run(
                plan=args.plan,
                paths=args.paths,
                analysis_workers=args.analysis_workers,
                process_existing=not args.no_process_existing,
                monitor=not args.no_monitor,
                watch=args.watch,
                patterns=args.patterns,
                ignore_patterns=args.ignore,
                recursive=not args.no_recursive,
                preset=args.preset,
                workflow_runner=None,
                center=args.center,
            )
        )
    except KeyboardInterrupt:
        pass
        