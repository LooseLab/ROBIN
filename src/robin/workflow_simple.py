"""Simple workflow management system for robin using threading."""

import os
import time
import threading
import queue
import itertools
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Any, Set

from watchdog.events import FileSystemEventHandler
from watchdog.observers import Observer
from tqdm import tqdm

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

    _RICH_AVAILABLE = True
except Exception:
    _RICH_AVAILABLE = False
from robin.logging_config import get_job_logger
from robin.analysis.target_analysis import (
    igv_bam_handler,
    snp_analysis_handler,
    target_bam_finalize_handler,
)

# ---------- Batch Configuration ----------
# See workflow_ray.py for the full description. Two timeouts are honoured:
#   - timeout_seconds:       used when no jobs of this type are in flight on
#                            a worker (idle). Keeps the first batch reactive.
#   - timeout_seconds_busy:  used when at least one job of this type is in
#                            flight, so the batcher accumulates more files
#                            into the next batch instead of dribbling out
#                            many tiny batches.
# Both can be globally overridden via env vars:
#   ROBIN_BATCH_TIMEOUT_IDLE_S, ROBIN_BATCH_TIMEOUT_BUSY_S
# Per-type override:
#   ROBIN_BATCH_TIMEOUT_BUSY_S_<TYPE>   e.g. ROBIN_BATCH_TIMEOUT_BUSY_S_CNV=45
BATCH_CONFIG: Dict[str, Dict[str, Any]] = {
    # Preprocessing should NOT be batched - each file needs individual sample ID extraction
    "preprocessing": {"max_batch_size": 1, "timeout_seconds": 0, "timeout_seconds_busy": 0},
    "bed_conversion": {"max_batch_size": 20, "timeout_seconds": 2, "timeout_seconds_busy": 30},
    "mgmt": {"max_batch_size": 20, "timeout_seconds": 2, "timeout_seconds_busy": 30},
    "cnv": {"max_batch_size": 50, "timeout_seconds": 2, "timeout_seconds_busy": 30},
    "target": {"max_batch_size": 20, "timeout_seconds": 2, "timeout_seconds_busy": 30},
    "fusion": {"max_batch_size": 20, "timeout_seconds": 2, "timeout_seconds_busy": 30},
    "sturgeon": {"max_batch_size": 20, "timeout_seconds": 2, "timeout_seconds_busy": 30},
    "nanodx": {"max_batch_size": 20, "timeout_seconds": 2, "timeout_seconds_busy": 30},
    "pannanodx": {"max_batch_size": 20, "timeout_seconds": 2, "timeout_seconds_busy": 30},
    "random_forest": {"max_batch_size": 20, "timeout_seconds": 2, "timeout_seconds_busy": 30},
    "igv_bam": {"max_batch_size": 20, "timeout_seconds": 2, "timeout_seconds_busy": 30},
    "snp_analysis": {"max_batch_size": 20, "timeout_seconds": 2, "timeout_seconds_busy": 30},
}


def _resolved_batch_timeouts(job_type: str) -> "tuple[float, float]":
    """Return (idle_timeout_s, busy_timeout_s) for a job_type, applying env overrides.

    Env precedence (most specific wins):
      1. ROBIN_BATCH_TIMEOUT_BUSY_S_<TYPE>
      2. ROBIN_BATCH_TIMEOUT_IDLE_S_<TYPE>
      3. ROBIN_BATCH_TIMEOUT_BUSY_S
      4. ROBIN_BATCH_TIMEOUT_IDLE_S
      5. BATCH_CONFIG defaults
    """
    cfg = BATCH_CONFIG.get(job_type, {})
    idle_default = float(cfg.get("timeout_seconds", 2) or 0)
    busy_default = float(cfg.get("timeout_seconds_busy", idle_default) or idle_default)

    def _f(name: str, fallback: float) -> float:
        try:
            v = os.environ.get(name)
            return float(v) if v is not None and v != "" else fallback
        except Exception:
            return fallback

    idle = _f("ROBIN_BATCH_TIMEOUT_IDLE_S", idle_default)
    busy = _f("ROBIN_BATCH_TIMEOUT_BUSY_S", busy_default)
    type_key = job_type.upper()
    idle = _f(f"ROBIN_BATCH_TIMEOUT_IDLE_S_{type_key}", idle)
    busy = _f(f"ROBIN_BATCH_TIMEOUT_BUSY_S_{type_key}", busy)
    if busy < idle:
        busy = idle
    return idle, busy


# === Shared Context for Each File ===
@dataclass
class WorkflowContext:
    """Context object that tracks metadata and results for each file being processed."""

    filepath: str
    metadata: dict = field(default_factory=dict)
    results: dict = field(default_factory=dict)
    history: list = field(default_factory=list)
    errors: list = field(default_factory=list)
    
    # NEW: Simple batch metadata
    batch_id: Optional[str] = None
    batch_index: Optional[int] = None

    def add_metadata(self, key: str, value: Any) -> None:
        """Add metadata to the context."""
        self.metadata[key] = value

    def add_result(self, job_type: str, result: Any) -> None:
        """Add a result from a job to the context."""
        self.results[job_type] = result
        self.history.append(job_type)

    def add_error(self, job_type: str, error: str) -> None:
        """Add an error to the context."""
        self.errors.append(
            {"job_type": job_type, "error": error, "timestamp": time.time()}
        )

    def get_summary(self) -> dict:
        """Get a summary of the workflow execution."""
        return {
            "filepath": self.filepath,
            "metadata": self.metadata,
            "results": self.results,
            "history": self.history,
            "errors": self.errors,
            "success": len(self.errors) == 0,
        }

    def get_sample_id(self) -> str:
        """Get the sample ID from metadata."""
        # First try to get from bam_metadata (for backward compatibility)
        bam_metadata = self.metadata.get("bam_metadata", {})
        sample_id = bam_metadata.get("sample_id", "unknown")
        
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

    def clear_all(self) -> None:
        """Clear all metadata/results/errors/history for this context.
        Use after all jobs for this file are finished to allow GC to reclaim memory.
        """
        self.metadata.clear()
        self.results.clear()
        self.history.clear()
        self.errors.clear()


# === Job Object ===
@dataclass
class Job:
    """Represents a single job in the workflow."""

    job_id: int
    job_type: str
    context: WorkflowContext
    origin: str  # "fast" or "slow"
    workflow: List[str]
    step: int = 0
    dependencies: Set[str] = field(default_factory=set)  # Job types this job depends on
    triggers: Set[str] = field(default_factory=set)  # Job types this job can trigger

    def next_job(self) -> Optional["Job"]:
        """Create the next job in the workflow if available."""
        if self.step + 1 < len(self.workflow):
            next_step = self.workflow[self.step + 1]
            # Ensure we have a string to work with
            if not isinstance(next_step, str):
                raise ValueError(
                    f"Expected string in workflow step, got {type(next_step)}: {next_step}"
                )

            if ":" not in next_step:
                raise ValueError(
                    f"Invalid workflow step format (missing ':'): {next_step}"
                )

            queue_type, next_type = next_step.split(":", 1)
            return Job(
                job_id=self.job_id,
                job_type=next_type,
                context=self.context,
                origin=queue_type,
                workflow=self.workflow,
                step=self.step + 1,
            )
        return None

    def get_sample_id(self) -> str:
        """Get the sample ID for this job."""
        return self.context.get_sample_id()

    def can_trigger_jobs(self) -> List["Job"]:
        """Create all jobs that can be triggered by this job's completion."""
        triggered_jobs = []

        # Define job dependencies and triggers
        job_dependencies = {
            # Jobs that depend on preprocessing
            "bed_conversion": {"preprocessing"},
            "mgmt": {"preprocessing"},
            "cnv": {"preprocessing"},
            "target": {"preprocessing"},
            "fusion": {"preprocessing"},
            # Jobs that depend on bed_conversion
            "sturgeon": {"bed_conversion"},
            "nanodx": {"bed_conversion"},
            "pannanodx": {"bed_conversion"},
            "random_forest": {"bed_conversion"},
        }

        # Get the original workflow plan to check what jobs are actually requested
        original_workflow = self.context.metadata.get("original_workflow", [])
        if original_workflow:
            # Extract job types from the workflow plan
            workflow_job_types = []
            for step in original_workflow:
                if isinstance(step, str) and ":" in step:
                    queue_type, job_type = step.split(":", 1)
                    workflow_job_types.append(job_type)
            original_workflow = workflow_job_types
        else:
            # Fallback: extract from current workflow
            original_workflow = []
            for step in self.workflow:
                if isinstance(step, str) and ":" in step:
                    queue_type, job_type = step.split(":", 1)
                    original_workflow.append(job_type)

        # Check if this job can trigger other jobs
        for job_type, dependencies in job_dependencies.items():
            if self.job_type in dependencies:
                # This job completion can trigger the dependent job
                # Check if all dependencies are met
                all_deps_met = True
                for dep in dependencies:
                    if dep not in self.context.history and dep != self.job_type:
                        all_deps_met = False
                        break

                # Only trigger if the job is in the original workflow plan
                if all_deps_met and job_type in original_workflow:
                    # Create the triggered job
                    queue_type = self._get_queue_type_for_job(job_type)
                    triggered_job = Job(
                        job_id=next(_job_id_counter),
                        job_type=job_type,
                        context=self.context,
                        origin=queue_type,
                        workflow=[f"{queue_type}:{job_type}"],  # Single step workflow
                        step=0,
                        dependencies=dependencies,
                        triggers=set(),  # No further triggers for now
                    )
                    triggered_jobs.append(triggered_job)

        return triggered_jobs

    def _get_queue_type_for_job(self, job_type: str) -> str:
        """Get the queue type for a given job type."""
        queue_mapping = {
            "preprocessing": "preprocessing",
            "bed_conversion": "bed_conversion",
            "mgmt": "mgmt",
            "cnv": "cnv",
            "target": "target",
            "fusion": "fusion",
            "sturgeon": "classification",
            "nanodx": "classification",
            "pannanodx": "classification",
            "random_forest": "slow",
        }
        return queue_mapping.get(job_type, "slow")


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


# === Sample Job Batcher ===
class SampleJobBatcher:
    def __init__(self, inflight_callback: Optional[Callable[[str], int]] = None):
        # sample_id -> job_type -> List[Job]
        self.pending_jobs: Dict[str, Dict[str, List[Job]]] = {}
        # sample_id -> job_type -> timestamp
        self.last_job_time: Dict[str, Dict[str, float]] = {}
        self.lock = threading.Lock()
        # Optional callback returning current inflight count for a job_type.
        # When > 0, the batcher uses the longer "busy" timeout so it
        # accumulates more files into one batch while the worker is occupied.
        self._inflight_callback: Optional[Callable[[str], int]] = inflight_callback

    def set_inflight_callback(self, cb: Optional[Callable[[str], int]]) -> None:
        """(Re)wire the inflight count callback after construction."""
        self._inflight_callback = cb

    def _inflight_for(self, job_type: str) -> int:
        if self._inflight_callback is None:
            return 0
        try:
            return int(self._inflight_callback(job_type) or 0)
        except Exception:
            return 0

    def _effective_timeout(self, job_type: str) -> float:
        """Pick idle vs busy timeout based on whether the type currently has
        any inflight jobs on a worker."""
        idle, busy = _resolved_batch_timeouts(job_type)
        return busy if self._inflight_for(job_type) > 0 else idle

    def add_job(self, job: Job) -> List[BatchedJob]:
        """Add a job and return any completed batches"""
        sample_id = job.get_sample_id()
        job_type = job.job_type
        
        with self.lock:
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
            return self._check_and_create_batches(sample_id, job_type)
    
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
    
    def check_timeouts(self) -> List[BatchedJob]:
        """Check for timed-out batches and return them"""
        current_time = time.time()
        timed_out_batches = []
        
        with self.lock:
            for sample_id in list(self.pending_jobs.keys()):
                for job_type in list(self.pending_jobs[sample_id].keys()):
                    jobs = self.pending_jobs[sample_id][job_type]
                    if not jobs:
                        continue
                    
                    timeout_seconds = self._effective_timeout(job_type)
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

    def force_flush_type(self, job_type: str) -> List[BatchedJob]:
        """Immediately flush whatever is pending for a given job_type across
        all samples, regardless of timeout. Called when a worker for that type
        becomes free so it has something to pick up instead of waiting out the
        busy timeout."""
        flushed: List[BatchedJob] = []
        with self.lock:
            for sample_id in list(self.pending_jobs.keys()):
                jobs = self.pending_jobs.get(sample_id, {}).get(job_type, [])
                if not jobs:
                    continue
                batch = self._create_batched_job(jobs, sample_id, job_type)
                flushed.append(batch)
                self.pending_jobs[sample_id][job_type] = []
                if (
                    sample_id in self.last_job_time
                    and job_type in self.last_job_time[sample_id]
                ):
                    del self.last_job_time[sample_id][job_type]
            self._cleanup_empty_entries()
        return flushed
    
    def _create_batched_job(self, jobs: List[Job], sample_id: str, job_type: str) -> BatchedJob:
        """Create a batched job from a list of individual jobs"""
        if not jobs:
            raise ValueError("Cannot create batched job from empty job list")
        
        # Validate all jobs have same sample_id and job_type
        for job in jobs:
            if job.get_sample_id() != sample_id:
                raise ValueError(f"Mixed sample IDs in batch: {sample_id} vs {job.get_sample_id()}")
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


# === Workflow Manager ===
class WorkflowManager:
    """Manages job queues and execution using threading with specialized workers."""

    def __init__(
        self,
        verbose: bool = False,
        analysis_workers: int = 1,
        use_separate_analysis_queues: bool = True,
        preprocessing_workers: int = 1,
        bed_conversion_workers: int = 1,
        enable_batching: bool = True,
    ):
        # Specialized queues for different task categories
        self.preprocessing_queue = queue.Queue()  # bam_preprocessing only
        self.bed_conversion_queue = queue.Queue()  # bed_conversion only

        # Analysis queue mode
        self.use_separate_analysis_queues = use_separate_analysis_queues

        if use_separate_analysis_queues:
            # Separate queues for each analysis job type (ensures only one of each type runs at a time)
            self.mgmt_queue = queue.Queue()  # MGMT analysis only
            self.cnv_queue = queue.Queue()  # CNV analysis only
            self.target_queue = queue.Queue()  # Target analysis only
            self.fusion_queue = queue.Queue()  # Fusion analysis only
            self.analysis_queue = None  # Not used in separate mode
        else:
            # Legacy single analysis queue (all analysis types share one queue)
            self.analysis_queue = queue.Queue()  # mgmt, cnv, target, fusion
            self.mgmt_queue = None
            self.cnv_queue = None
            self.target_queue = None
            self.fusion_queue = None

        self.classification_queue = queue.Queue()  # sturgeon, nanodx, pannanodx
        self.slow_queue = queue.Queue()  # slow jobs (legacy)

        self.completed_jobs = []
        self.failed_jobs = []

        # Analysis workers per queue (each analysis type gets its own worker)
        self.analysis_workers_count = analysis_workers
        # Configurable workers for preprocessing and bed conversion (defaults to 1)
        self.preprocessing_workers_count = (
            preprocessing_workers
            if isinstance(preprocessing_workers, int) and preprocessing_workers > 0
            else 1
        )
        self.bed_conversion_workers_count = (
            bed_conversion_workers
            if isinstance(bed_conversion_workers, int) and bed_conversion_workers > 0
            else 1
        )
        self.classification_workers_count = 1
        self.slow_workers_count = 1

        # Specialized job handlers for each queue
        self.job_handlers_preprocessing: Dict[str, Callable[[Job], None]] = {}
        self.job_handlers_bed_conversion: Dict[str, Callable[[Job], None]] = {}

        # Analysis job handlers (mode-dependent)
        if use_separate_analysis_queues:
            # Separate handlers for each analysis job type
            self.job_handlers_mgmt: Dict[str, Callable[[Job], None]] = {}
            self.job_handlers_cnv: Dict[str, Callable[[Job], None]] = {}
            self.job_handlers_target: Dict[str, Callable[[Job], None]] = {}
            self.job_handlers_fusion: Dict[str, Callable[[Job], None]] = {}
            self.job_handlers_analysis = None  # Not used in separate mode
        else:
            # Legacy single analysis handler
            self.job_handlers_analysis: Dict[str, Callable[[Job], None]] = {}
            self.job_handlers_mgmt = None
            self.job_handlers_cnv = None
            self.job_handlers_target = None
            self.job_handlers_fusion = None

        self.job_handlers_classification: Dict[str, Callable[[Job], None]] = {}
        self.job_handlers_slow: Dict[str, Callable[[Job], None]] = {}

        self.running = True
        self.verbose = verbose

        # Job deduplication tracking
        # Jobs that should be deduplicated by sample ID (e.g., sturgeon analysis)
        self.deduplicate_job_types: Set[str] = {
            "sturgeon",
            "nanodx",
            "pannanodx",
            "random_forest",
        }
        # Track running and pending jobs by sample ID for deduplication
        # Allow max 2 jobs per sample: 1 running + 1 pending
        self.running_jobs_by_sample: Dict[str, Dict[str, int]] = (
            {}
        )  # sample_id -> {job_type: count}
        self.pending_jobs_by_sample: Dict[str, Dict[str, int]] = (
            {}
        )  # sample_id -> {job_type: count}

        self.job_tracking_lock = threading.Lock()

        # Track worker threads for graceful shutdown
        self.preprocessing_workers = []
        self.bed_conversion_workers = []

        # Analysis workers (mode-dependent)
        if use_separate_analysis_queues:
            # Separate workers for each analysis job type
            self.mgmt_workers = []
            self.cnv_workers = []
            self.target_workers = []
            self.fusion_workers = []
            self.analysis_workers = []  # Not used in separate mode
        else:
            # Legacy single analysis workers
            self.analysis_workers = []
            self.mgmt_workers = []
            self.cnv_workers = []
            self.target_workers = []
            self.fusion_workers = []

        self.classification_workers = []
        self.slow_workers = []

        # Progress tracking
        self.total_jobs_enqueued = 0
        self.total_jobs_skipped = 0  # Track jobs skipped due to deduplication
        self.active_jobs = {}  # job_id -> job_info
        self.job_start_times = {}  # job_id -> start_time
        self.completed_jobs_by_type = (
            {}
        )  # job_type -> count (successful completions only)
        self.failed_jobs_by_type = {}  # job_type -> count (failed jobs only)
        self.progress_lock = threading.Lock()  # Lock for progress tracking

        # Per-sample tracking for GUI/monitoring
        # sample_id -> {
        #   'sample_id', 'active_jobs', 'total_jobs', 'completed_jobs', 'failed_jobs', 'job_types' (set), 'last_seen'
        # }
        self.samples_by_id: Dict[str, Dict[str, Any]] = {}
        
        # Batching support. The batcher uses adaptive timeouts driven by
        # the live per-type inflight count from self.active_jobs. When a
        # type is idle, fire batches quickly; when busy, accumulate longer.
        self.enable_batching = enable_batching
        self.sample_batcher = (
            SampleJobBatcher(
                inflight_callback=self._inflight_count_for_type
            )
            if enable_batching
            else None
        )

        # Coalesce stats (exposed via get_stats() for debugging)
        self.coalesce_events: int = 0
        self.coalesce_jobs_absorbed: int = 0
        self.coalesce_contexts_merged: int = 0

        # Define job type mappings to queues
        if use_separate_analysis_queues:
            self.job_queue_mapping = {
                # Preprocessing queue
                "preprocessing": "preprocessing",
                # Bed conversion queue
                "bed_conversion": "bed_conversion",
                # Separate analysis queues (ensures only one of each type runs at a time)
                "mgmt": "mgmt",
                "cnv": "cnv",
                "target": "target",
                "fusion": "fusion",
                "test": "mgmt",  # For testing purposes - map to mgmt queue
                "long": "cnv",  # For testing purposes - map to cnv queue
                "quick": "target",  # For testing purposes - map to target queue
                # IGV BAM build goes to slow queue to avoid contention with analyses
                "igv_bam": "slow",
                "target_bam_finalize": "slow",
                # Classification queue
                "sturgeon": "classification",
                "nanodx": "classification",
                "pannanodx": "classification",
                # Slow queue
                "random_forest": "slow",
                "sleep": "slow",
                "echo": "slow",
            }
        else:
            self.job_queue_mapping = {
                # Preprocessing queue
                "preprocessing": "preprocessing",
                # Bed conversion queue
                "bed_conversion": "bed_conversion",
                # Legacy single analysis queue
                "mgmt": "analysis",
                "cnv": "analysis",
                "target": "analysis",
                "fusion": "analysis",
                "test": "analysis",  # For testing purposes
                "long": "analysis",  # For testing purposes
                "quick": "analysis",  # For testing purposes
                "igv_bam": "slow",
                "target_bam_finalize": "slow",
                # Classification queue
                "sturgeon": "classification",
                "nanodx": "classification",
                "pannanodx": "classification",
                # Slow queue
                "random_forest": "slow",
                "sleep": "slow",
                "echo": "slow",
            }

    def register_handler(
        self, queue_type: str, job_type: str, handler: Callable[[Job], None]
    ) -> None:
        """Register a handler for a specific job type and queue."""
        if queue_type == "preprocessing":
            self.job_handlers_preprocessing[job_type] = handler
        elif queue_type == "bed_conversion":
            self.job_handlers_bed_conversion[job_type] = handler
        elif queue_type == "analysis":
            if not self.use_separate_analysis_queues:
                self.job_handlers_analysis[job_type] = handler
            else:
                raise ValueError(
                    "Cannot register 'analysis' queue handler in separate analysis queues mode"
                )
        elif queue_type == "mgmt":
            if self.use_separate_analysis_queues:
                self.job_handlers_mgmt[job_type] = handler
            else:
                raise ValueError(
                    "Cannot register 'mgmt' queue handler in legacy analysis queue mode"
                )
        elif queue_type == "cnv":
            if self.use_separate_analysis_queues:
                self.job_handlers_cnv[job_type] = handler
            else:
                raise ValueError(
                    "Cannot register 'cnv' queue handler in legacy analysis queue mode"
                )
        elif queue_type == "target":
            if self.use_separate_analysis_queues:
                self.job_handlers_target[job_type] = handler
            else:
                raise ValueError(
                    "Cannot register 'target' queue handler in legacy analysis queue mode"
                )
        elif queue_type == "fusion":
            if self.use_separate_analysis_queues:
                self.job_handlers_fusion[job_type] = handler
            else:
                raise ValueError(
                    "Cannot register 'fusion' queue handler in legacy analysis queue mode"
                )
        elif queue_type == "classification":
            self.job_handlers_classification[job_type] = handler
        elif queue_type == "slow":
            self.job_handlers_slow[job_type] = handler
        else:
            raise ValueError(f"Invalid queue type: {queue_type}")

    def add_deduplication_job_type(self, job_type: str) -> None:
        """Add a job type to the deduplication list."""
        self.deduplicate_job_types.add(job_type)

    def remove_deduplication_job_type(self, job_type: str) -> None:
        """Remove a job type from the deduplication list."""
        self.deduplicate_job_types.discard(job_type)

    def _can_enqueue_job_for_sample(self, job_type: str, sample_id: str) -> bool:
        """Check if a job can be enqueued for the sample (max 2: 1 running + 1 pending)."""
        with self.job_tracking_lock:
            running_count = 0
            pending_count = 0

            if (
                sample_id in self.running_jobs_by_sample
                and job_type in self.running_jobs_by_sample[sample_id]
            ):
                running_count = self.running_jobs_by_sample[sample_id][job_type]

            if (
                sample_id in self.pending_jobs_by_sample
                and job_type in self.pending_jobs_by_sample[sample_id]
            ):
                pending_count = self.pending_jobs_by_sample[sample_id][job_type]

            # Allow if we have less than 2 total jobs (running + pending)
            # AND we don't exceed 1 pending job
            return (running_count + pending_count) < 2 and pending_count < 1

    def _mark_job_pending_for_sample(self, job_type: str, sample_id: str) -> None:
        """Mark a job as pending for the sample."""
        with self.job_tracking_lock:
            if sample_id not in self.pending_jobs_by_sample:
                self.pending_jobs_by_sample[sample_id] = {}
            if job_type not in self.pending_jobs_by_sample[sample_id]:
                self.pending_jobs_by_sample[sample_id][job_type] = 0
            self.pending_jobs_by_sample[sample_id][job_type] += 1

    def _mark_job_running_for_sample(self, job_type: str, sample_id: str) -> None:
        """Mark a job as running for the sample (move from pending to running)."""
        with self.job_tracking_lock:
            # Remove from pending
            if (
                sample_id in self.pending_jobs_by_sample
                and job_type in self.pending_jobs_by_sample[sample_id]
            ):
                self.pending_jobs_by_sample[sample_id][job_type] -= 1
                if self.pending_jobs_by_sample[sample_id][job_type] <= 0:
                    del self.pending_jobs_by_sample[sample_id][job_type]
                if not self.pending_jobs_by_sample[sample_id]:
                    del self.pending_jobs_by_sample[sample_id]

            # Add to running
            if sample_id not in self.running_jobs_by_sample:
                self.running_jobs_by_sample[sample_id] = {}
            if job_type not in self.running_jobs_by_sample[sample_id]:
                self.running_jobs_by_sample[sample_id][job_type] = 0
            self.running_jobs_by_sample[sample_id][job_type] += 1

    def _unmark_job_for_sample(self, job_type: str, sample_id: str) -> None:
        """Unmark a job as running for the sample (job completed)."""
        with self.job_tracking_lock:
            # Remove from running
            if (
                sample_id in self.running_jobs_by_sample
                and job_type in self.running_jobs_by_sample[sample_id]
            ):
                self.running_jobs_by_sample[sample_id][job_type] -= 1
                if self.running_jobs_by_sample[sample_id][job_type] <= 0:
                    del self.running_jobs_by_sample[sample_id][job_type]
                if not self.running_jobs_by_sample[sample_id]:
                    del self.running_jobs_by_sample[sample_id]

    def enqueue_jobs(self, jobs: List[Job]) -> None:
        """Add jobs to their respective queues with bounded deduplication (max 2: 1 running + 1 pending)."""
        # Process jobs through batcher if batching is enabled
        if self.enable_batching and self.sample_batcher:
            # Separate preprocessing jobs from other jobs
            preprocessing_jobs = [job for job in jobs if job.job_type == "preprocessing"]
            other_jobs = [job for job in jobs if job.job_type != "preprocessing"]
            
            # Submit preprocessing jobs directly (no batching)
            if preprocessing_jobs:
                self._enqueue_jobs_internal(preprocessing_jobs)
            
            # Process other jobs through batcher
            for job in other_jobs:
                batches = self.sample_batcher.add_job(job)
                if batches:
                    # Convert and enqueue batches
                    regular_jobs = self._convert_batched_jobs(batches)
                    self._enqueue_jobs_internal(regular_jobs)
        else:
            self._enqueue_jobs_internal(jobs)
    
    def _enqueue_jobs_internal(self, jobs: List[Job]) -> None:
        """Internal job enqueue logic"""
        jobs_to_enqueue = []

        for job in jobs:
            sample_id = job.get_sample_id()

            # Check if this job type should be deduplicated
            if job.job_type in self.deduplicate_job_types:
                if not self._can_enqueue_job_for_sample(job.job_type, sample_id):
                    # Already have max jobs (1 running + 1 pending) for this sample, skip it
                    logger = get_job_logger(
                        str(job.job_id), job.job_type, job.context.filepath
                    )
                    logger.info(
                        f"Skipping {job.job_type} job for sample {sample_id} (max jobs reached: 1 running + 1 pending)"
                    )
                    self.total_jobs_skipped += 1  # Track skipped jobs
                    continue
                else:
                    # Mark this job as pending for the sample
                    self._mark_job_pending_for_sample(job.job_type, sample_id)
                    logger = get_job_logger(
                        str(job.job_id), job.job_type, job.context.filepath
                    )
                    logger.debug(
                        f"Marked {job.job_type} job as pending for sample {sample_id}"
                    )

            jobs_to_enqueue.append(job)

        # Enqueue the filtered jobs
        for job in jobs_to_enqueue:
            # Track total jobs enqueued
            self.total_jobs_enqueued += 1

            # Determine which queue to use based on job type
            queue_type = self.job_queue_mapping.get(job.job_type, "slow")

            # Log job enqueuing for debugging
            logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
            logger.info(
                f"Enqueuing job {job.job_id} ({job.job_type}) to {queue_type} queue"
            )

            if queue_type == "preprocessing":
                self.preprocessing_queue.put(job)
            elif queue_type == "bed_conversion":
                self.bed_conversion_queue.put(job)
            elif queue_type == "analysis":
                if not self.use_separate_analysis_queues:
                    self.analysis_queue.put(job)
                else:
                    raise ValueError(
                        "Cannot enqueue to 'analysis' queue in separate analysis queues mode"
                    )
            elif queue_type == "mgmt":
                if self.use_separate_analysis_queues:
                    self.mgmt_queue.put(job)
                else:
                    raise ValueError(
                        "Cannot enqueue to 'mgmt' queue in legacy analysis queue mode"
                    )
            elif queue_type == "cnv":
                if self.use_separate_analysis_queues:
                    self.cnv_queue.put(job)
                else:
                    raise ValueError(
                        "Cannot enqueue to 'cnv' queue in legacy analysis queue mode"
                    )
            elif queue_type == "target":
                if self.use_separate_analysis_queues:
                    self.target_queue.put(job)
                else:
                    raise ValueError(
                        "Cannot enqueue to 'target' queue in legacy analysis queue mode"
                    )
            elif queue_type == "fusion":
                if self.use_separate_analysis_queues:
                    self.fusion_queue.put(job)
                else:
                    raise ValueError(
                        "Cannot enqueue to 'fusion' queue in legacy analysis queue mode"
                    )
            elif queue_type == "classification":
                self.classification_queue.put(job)
            elif queue_type == "slow":
                self.slow_queue.put(job)
            else:
                raise ValueError(
                    f"Invalid queue type for job {job.job_type}: {queue_type}"
                )

    def worker(
        self, name: str, job_queue: queue.Queue, handlers: Dict[str, Callable]
    ) -> None:
        """Worker that processes jobs from a queue."""
        while self.running:
            try:
                job = job_queue.get(timeout=1.0)

                # Try to absorb compatible waiting siblings into this job's
                # batch before processing. Mops up small-batch fragmentation
                # caused by the per-(sample, type) batcher firing on its
                # busy/idle timeout while the worker was occupied. No-op for
                # non-batched jobs and for jobs tagged force_individual_batch.
                if self.enable_batching:
                    try:
                        sid_for_coalesce = (
                            job.get_sample_id()
                            if hasattr(job, "get_sample_id")
                            else "unknown"
                        )
                        if sid_for_coalesce and sid_for_coalesce != "unknown":
                            self._coalesce_from_queue(
                                job,
                                job_queue,
                                sid_for_coalesce,
                                job.job_type,
                            )
                    except Exception:
                        pass

                # Track job start time and active status (this includes queue waiting time)
                self.job_start_times[job.job_id] = time.time()
                self.active_jobs[job.job_id] = {
                    "job_type": job.job_type,
                    "filepath": job.context.filepath,
                    "sample_id": (
                        job.get_sample_id()
                        if hasattr(job, "get_sample_id")
                        else "unknown"
                    ),
                    "worker": name,
                    "start_time": time.time(),
                    "processing_start_time": None,  # Will be set when actual processing starts
                }

                # Get job-specific logger
                logger = get_job_logger(
                    str(job.job_id), job.job_type, job.context.filepath
                )
                logger.info(
                    f"Starting job {job.job_id} ({job.job_type}) for {job.context.filepath}"
                )

                # Mark job as running for deduplicated job types
                if job.job_type in self.deduplicate_job_types:
                    sample_id = job.get_sample_id()
                    self._mark_job_running_for_sample(job.job_type, sample_id)
                    logger.debug(
                        f"Marked {job.job_type} job as running for sample {sample_id}"
                    )

                # Per-sample start tracking (when sample_id is known)
                sid_start = job.get_sample_id()
                if sid_start and sid_start != "unknown":
                    self._on_sample_job_started(sid_start, job.job_type)

                try:
                    # Set processing start time right before calling the handler (excludes queue waiting time)
                    self.active_jobs[job.job_id]["processing_start_time"] = time.time()

                    handlers[job.job_type](job)

                    # Track this job as completed regardless of whether it has a next job
                    logger.info(
                        f"Job {job.job_id} ({job.job_type}) completed successfully"
                    )
                    self.completed_jobs.append(job.job_id)

                    # Track completion by job type
                    with self.progress_lock:
                        if job.job_type not in self.completed_jobs_by_type:
                            self.completed_jobs_by_type[job.job_type] = 0
                        self.completed_jobs_by_type[job.job_type] += 1

                    # Check if this job can trigger other jobs (parallel execution)
                    triggered_jobs = job.can_trigger_jobs()
                    if triggered_jobs:
                        logger.info(
                            f"Job {job.job_type} completed, triggering {len(triggered_jobs)} parallel jobs: {[j.job_type for j in triggered_jobs]}"
                        )
                        
                        # For CNV jobs, use batching if enabled
                        if self.enable_batching and self.sample_batcher:
                            # Check if any of the triggered jobs are CNV jobs
                            cnv_jobs = [j for j in triggered_jobs if j.job_type == "cnv"]
                            other_jobs = [j for j in triggered_jobs if j.job_type != "cnv"]
                            
                            # Submit non-CNV jobs immediately
                            if other_jobs:
                                self.enqueue_jobs(other_jobs)
                            
                            # For CNV jobs, add them to the batcher
                            for cnv_job in cnv_jobs:
                                batches = self.sample_batcher.add_job(cnv_job)
                                if batches:
                                    # Convert BatchedJob to regular Job for processing
                                    regular_jobs = self._convert_batched_jobs(batches)
                                    self._enqueue_jobs_internal(regular_jobs)
                        else:
                            # No batching, submit all jobs normally
                            self.enqueue_jobs(triggered_jobs)
                    else:
                        # Fallback to linear workflow if no parallel triggers
                        next_job = job.next_job()
                        if next_job:
                            logger.info(f"Creating next job: {next_job.job_type}")
                            self.enqueue_jobs([next_job])
                        else:
                            logger.info(
                                f"No more jobs in workflow for {job.context.filepath}"
                            )
                            logger.debug(f"Final results: {job.context.results}")
                            # Context clearing temporarily disabled to avoid interfering with downstream queues

                except Exception as e:
                    error_msg = f"Job {job.job_type} failed: {str(e)}"
                    job.context.add_error(job.job_type, error_msg)
                    self.failed_jobs.append(job.job_id)

                    # Track failed jobs by type
                    with self.progress_lock:
                        if job.job_type not in self.failed_jobs_by_type:
                            self.failed_jobs_by_type[job.job_type] = 0
                        self.failed_jobs_by_type[job.job_type] += 1

                    logger.error(f"Job {job.job_id} failed: {error_msg}")
                finally:
                    # Per-sample finish tracking (prefer sample_id after handler updated context)
                    sid_finish = job.get_sample_id()
                    if sid_finish and sid_finish != "unknown":
                        self._on_sample_job_finished(
                            sid_finish, job.job_type, job.job_id in self.completed_jobs
                        )

                    # Remove from active jobs (this lowers our inflight count
                    # for this type, so the batcher will start using the idle
                    # timeout for it again).
                    if job.job_id in self.active_jobs:
                        del self.active_jobs[job.job_id]
                    if job.job_id in self.job_start_times:
                        del self.job_start_times[job.job_id]

                    # Always unmark the job as running when it's done (success or failure)
                    if job.job_type in self.deduplicate_job_types:
                        sample_id = job.get_sample_id()
                        self._unmark_job_for_sample(job.job_type, sample_id)
                        logger.debug(
                            f"Unmarked {job.job_type} job as running for sample {sample_id}"
                        )

                    # A worker for this type just freed up. Force-flush any
                    # batches the batcher is sitting on so this worker (or
                    # another idle worker) has something to pick up
                    # immediately, without waiting out the busy timeout.
                    try:
                        self._flush_batcher_for_type(job.job_type)
                    except Exception:
                        pass

                job_queue.task_done()

            except queue.Empty:
                continue

    def run(self) -> None:
        """Start the workflow manager with specialized workers for different task categories."""
        # Create preprocessing workers (bam_preprocessing only)
        for i in range(self.preprocessing_workers_count):
            preprocessing_worker = threading.Thread(
                target=self.worker,
                args=(
                    f"PreprocessingWorker-{i+1}",
                    self.preprocessing_queue,
                    self.job_handlers_preprocessing,
                ),
            )
            preprocessing_worker.daemon = True
            self.preprocessing_workers.append(preprocessing_worker)
            preprocessing_worker.start()

        # Create bed conversion workers (bed_conversion only)
        for i in range(self.bed_conversion_workers_count):
            bed_conversion_worker = threading.Thread(
                target=self.worker,
                args=(
                    f"BedConversionWorker-{i+1}",
                    self.bed_conversion_queue,
                    self.job_handlers_bed_conversion,
                ),
            )
            bed_conversion_worker.daemon = True
            self.bed_conversion_workers.append(bed_conversion_worker)
            bed_conversion_worker.start()

        # Create analysis workers based on mode
        if self.use_separate_analysis_queues:
            # Create separate workers for each analysis job type (ensures only one of each type runs at a time)
            for i in range(self.analysis_workers_count):
                # MGMT worker
                mgmt_worker = threading.Thread(
                    target=self.worker,
                    args=(f"MGMTWorker-{i+1}", self.mgmt_queue, self.job_handlers_mgmt),
                )
                mgmt_worker.daemon = True
                self.mgmt_workers.append(mgmt_worker)
                mgmt_worker.start()

                # CNV worker
                cnv_worker = threading.Thread(
                    target=self.worker,
                    args=(f"CNVWorker-{i+1}", self.cnv_queue, self.job_handlers_cnv),
                )
                cnv_worker.daemon = True
                self.cnv_workers.append(cnv_worker)
                cnv_worker.start()

                # Target worker
                target_worker = threading.Thread(
                    target=self.worker,
                    args=(
                        f"TargetWorker-{i+1}",
                        self.target_queue,
                        self.job_handlers_target,
                    ),
                )
                target_worker.daemon = True
                self.target_workers.append(target_worker)
                target_worker.start()

                # Fusion worker
                fusion_worker = threading.Thread(
                    target=self.worker,
                    args=(
                        f"FusionWorker-{i+1}",
                        self.fusion_queue,
                        self.job_handlers_fusion,
                    ),
                )
                fusion_worker.daemon = True
                self.fusion_workers.append(fusion_worker)
                fusion_worker.start()
        else:
            # Create legacy single analysis workers
            for i in range(self.analysis_workers_count):
                analysis_worker = threading.Thread(
                    target=self.worker,
                    args=(
                        f"AnalysisWorker-{i+1}",
                        self.analysis_queue,
                        self.job_handlers_analysis,
                    ),
                )
                analysis_worker.daemon = True
                self.analysis_workers.append(analysis_worker)
                analysis_worker.start()

        # Create classification workers (sturgeon, nanodx, pannanodx)
        for i in range(self.classification_workers_count):
            classification_worker = threading.Thread(
                target=self.worker,
                args=(
                    f"ClassificationWorker-{i+1}",
                    self.classification_queue,
                    self.job_handlers_classification,
                ),
            )
            classification_worker.daemon = True
            self.classification_workers.append(classification_worker)
            classification_worker.start()

        # Create slow workers (legacy)
        for i in range(self.slow_workers_count):
            slow_worker = threading.Thread(
                target=self.worker,
                args=(f"SlowWorker-{i+1}", self.slow_queue, self.job_handlers_slow),
            )
            slow_worker.daemon = True
            self.slow_workers.append(slow_worker)
            slow_worker.start()
        
        # Start batch timeout thread if batching is enabled
        if self.enable_batching and self.sample_batcher:
            self.batch_timeout_thread = threading.Thread(
                target=self._batch_timeout_loop, daemon=True
            )
            self.batch_timeout_thread.start()

    def _try_clear_context(self, context: "WorkflowContext") -> None:
        """Clear a context only when all jobs in the original workflow have produced results.
        This prevents premature clearing that would block downstream triggers.
        """
        original_workflow = context.metadata.get("original_workflow", [])
        required_job_types: Set[str] = set()
        for step in original_workflow:
            if isinstance(step, str) and ":" in step:
                _, job_type = step.split(":", 1)
                required_job_types.add(job_type)
        if not required_job_types:
            return
        completed_job_types = set(context.history)
        # Only clear when every planned job type has at least one result entry
        if required_job_types.issubset(completed_job_types):
            context.clear_all()

        # Wait for all workers to finish with periodic checks for shutdown
        all_workers = (
            self.preprocessing_workers
            + self.bed_conversion_workers
            + self.mgmt_workers
            + self.cnv_workers
            + self.target_workers
            + self.fusion_workers
            + self.classification_workers
            + self.slow_workers
        )

        try:
            while self.running:
                # Check if any workers are still alive
                alive_workers = [w for w in all_workers if w.is_alive()]
                if not alive_workers:
                    break

                # Sleep briefly to allow for interrupt handling
                time.sleep(0.1)

        except KeyboardInterrupt:
            print("\n[SHUTDOWN] Interrupted by user (Ctrl-C)")
            print("[SHUTDOWN] Initiating graceful shutdown...")
            # Signal workers to stop
            self.running = False
            print("[SHUTDOWN] Signalled workers to stop")
            raise

    def stop(self, timeout: float = 30.0) -> bool:
        """
        Stop the workflow manager and clean up running threads.

        Args:
            timeout: Maximum time to wait for workers to finish (seconds)

        Returns:
            True if all workers stopped gracefully, False if timeout occurred
        """
        if not self.running:
            return True

        if self.verbose:
            print(f"[WorkflowManager] Stopping workflow manager (timeout: {timeout}s)")

        # Signal workers to stop
        self.running = False

        # Wait for workers to finish with timeout
        if self.use_separate_analysis_queues:
            all_workers = (
                self.preprocessing_workers
                + self.bed_conversion_workers
                + self.mgmt_workers
                + self.cnv_workers
                + self.target_workers
                + self.fusion_workers
                + self.classification_workers
                + self.slow_workers
            )
        else:
            all_workers = (
                self.preprocessing_workers
                + self.bed_conversion_workers
                + self.analysis_workers
                + self.classification_workers
                + self.slow_workers
            )

        # Wait for each worker with timeout
        for worker in all_workers:
            if worker.is_alive():
                worker.join(timeout=timeout)
                if worker.is_alive():
                    if self.verbose:
                        print(
                            f"[WorkflowManager] Warning: Worker {worker.name} did not stop within timeout"
                        )
                    return False

        if self.verbose:
            print("[WorkflowManager] All workers stopped successfully")

        return True

    def shutdown(self, timeout: float = 30.0) -> bool:
        """
        Alias for stop() method for backward compatibility.

        Args:
            timeout: Maximum time to wait for workers to finish (seconds)

        Returns:
            True if all workers stopped gracefully, False if timeout occurred
        """
        return self.stop(timeout)
    
    def _batch_timeout_loop(self):
        """Periodically check for timed-out batches"""
        while self.running:
            try:
                # Check for timed-out batches
                timed_out_batches = self.sample_batcher.check_timeouts()
                
                # Enqueue timed-out batches
                if timed_out_batches:
                    regular_jobs = self._convert_batched_jobs(timed_out_batches)
                    self._enqueue_jobs_internal(regular_jobs)
                
                time.sleep(1.0)  # Check every second
                
            except Exception:
                time.sleep(1.0)

    # --- Adaptive timeout / coalesce helpers ---
    def _inflight_count_for_type(self, job_type: str) -> int:
        """Return the number of jobs of ``job_type`` currently being processed
        on a worker. Used by SampleJobBatcher to switch between idle and busy
        timeouts."""
        try:
            return sum(
                1 for info in self.active_jobs.values()
                if info.get("job_type") == job_type
            )
        except Exception:
            return 0

    def _flush_batcher_for_type(self, job_type: str) -> None:
        """Force-flush whatever the batcher is holding for ``job_type`` and
        push the resulting batches into the relevant queue. Called when a
        worker for that type frees up so it has something to pick up
        immediately without waiting out the busy timeout."""
        if not (self.enable_batching and self.sample_batcher):
            return
        try:
            flushed = self.sample_batcher.force_flush_type(job_type)
        except Exception:
            return
        if not flushed:
            return
        regular_jobs = self._convert_batched_jobs(flushed)
        try:
            self._enqueue_jobs_internal(regular_jobs)
        except Exception:
            pass

    @staticmethod
    def _job_is_coalescable_batch(job: Optional[Job]) -> bool:
        """A job is mergeable iff it carries a _batched_job and none of its
        contexts opted out via force_individual_batch."""
        if job is None:
            return False
        try:
            bjob = job.context.metadata.get("_batched_job")
        except Exception:
            return False
        if bjob is None:
            return False
        try:
            for ctx in bjob.contexts:
                if (ctx.metadata or {}).get("force_individual_batch"):
                    return False
        except Exception:
            return False
        return True

    def _coalesce_from_queue(
        self,
        job: Job,
        job_queue: queue.Queue,
        sample_id: str,
        job_type: str,
    ) -> int:
        """Absorb compatible siblings from ``job_queue`` into ``job``'s batch.

        Scans the queue for other Jobs with the same ``sample_id`` and
        ``job_type`` that also carry a _batched_job and aren't tagged
        force_individual_batch, and merges their contexts into this job's
        _batched_job up to ``max_batch_size``. Absorbed entries are removed
        from the queue. Returns the number of waiting Jobs absorbed.
        """
        if not self._job_is_coalescable_batch(job):
            return 0
        if not sample_id or sample_id == "unknown":
            return 0

        bjob: BatchedJob = job.context.metadata.get("_batched_job")
        cfg = BATCH_CONFIG.get(job_type, {"max_batch_size": 20})
        try:
            max_size = int(cfg.get("max_batch_size", 20))
        except Exception:
            max_size = 20
        if max_size <= 1:
            return 0
        room = max_size - len(bjob.contexts)
        if room <= 0:
            return 0

        absorbed = 0
        merged_contexts = 0

        # Manipulate the underlying deque under the queue's mutex. We adjust
        # unfinished_tasks for each removed entry so q.task_done() / q.join()
        # accounting stays balanced.
        with job_queue.mutex:
            dq = job_queue.queue
            i = 0
            while i < len(dq) and room > 0:
                entry = dq[i]
                if not isinstance(entry, Job) or entry.job_type != job_type:
                    i += 1
                    continue
                try:
                    entry_sid = entry.get_sample_id()
                except Exception:
                    i += 1
                    continue
                if entry_sid != sample_id:
                    i += 1
                    continue
                if not self._job_is_coalescable_batch(entry):
                    i += 1
                    continue
                other_bjob: BatchedJob = entry.context.metadata.get(
                    "_batched_job"
                )
                take = min(room, len(other_bjob.contexts))
                if take <= 0:
                    i += 1
                    continue
                bjob.contexts.extend(other_bjob.contexts[:take])
                room -= take
                merged_contexts += take
                if take >= len(other_bjob.contexts):
                    # Whole entry absorbed; remove from the queue and balance
                    # unfinished_tasks (which was incremented by put()).
                    del dq[i]
                    try:
                        if job_queue.unfinished_tasks > 0:
                            job_queue.unfinished_tasks -= 1
                            if (
                                job_queue.unfinished_tasks == 0
                                and hasattr(job_queue, "all_tasks_done")
                            ):
                                job_queue.all_tasks_done.notify_all()
                    except Exception:
                        pass
                    absorbed += 1
                else:
                    other_bjob.contexts = other_bjob.contexts[take:]
                    i += 1
                    # Loop will exit on room == 0.

        if absorbed > 0:
            try:
                bjob.batch_id = (
                    f"{sample_id}_{job_type}_coalesced_{int(time.time() * 1000)}"
                )
            except Exception:
                pass
            try:
                self.coalesce_events += 1
                self.coalesce_jobs_absorbed += absorbed
                self.coalesce_contexts_merged += merged_contexts
            except Exception:
                pass
        return absorbed
    
    def _convert_batched_jobs(self, batched_jobs: List[BatchedJob]) -> List[Job]:
        """Convert BatchedJob objects to regular Job objects for processing"""
        regular_jobs = []
        for batched_job in batched_jobs:
            # Create a regular job that will handle the batch
            regular_job = Job(
                job_id=batched_job.job_id,
                job_type=batched_job.job_type,
                context=batched_job.contexts[0],  # Use first context as primary
                origin=batched_job.origin,
                workflow=batched_job.workflow,
                step=batched_job.step
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

    def is_running(self) -> bool:
        """Check if the workflow manager is currently running."""
        return self.running

    def get_stats(self) -> dict:
        """Get statistics about the workflow execution."""
        # Get queue sizes
        preprocessing_queue_size = self.preprocessing_queue.qsize()
        bed_conversion_queue_size = self.bed_conversion_queue.qsize()

        if self.use_separate_analysis_queues:
            mgmt_queue_size = self.mgmt_queue.qsize()
            cnv_queue_size = self.cnv_queue.qsize()
            target_queue_size = self.target_queue.qsize()
            fusion_queue_size = self.fusion_queue.qsize()
            analysis_queue_size = 0  # Not used in separate mode
        else:
            analysis_queue_size = self.analysis_queue.qsize()
            mgmt_queue_size = 0
            cnv_queue_size = 0
            target_queue_size = 0
            fusion_queue_size = 0

        classification_queue_size = self.classification_queue.qsize()
        slow_queue_size = self.slow_queue.qsize()

        # Get active jobs by worker type
        active_by_worker = {}
        for job_info in self.active_jobs.values():
            worker = job_info["worker"]
            if worker not in active_by_worker:
                active_by_worker[worker] = []

            # Use processing_start_time if available (excludes queue waiting time), otherwise fall back to start_time
            if job_info.get("processing_start_time") is not None:
                duration = time.time() - job_info["processing_start_time"]
            else:
                duration = time.time() - job_info["start_time"]

            active_by_worker[worker].append(
                {
                    "job_type": job_info["job_type"],
                    "filepath": job_info["filepath"],
                    "sample_id": job_info.get("sample_id", "unknown"),
                    "duration": duration,
                }
            )

        # Calculate total expected jobs based on completed + failed + active + queued
        # This gives us the current total of all jobs that have been processed or are in progress
        if self.use_separate_analysis_queues:
            total_expected = (
                len(self.completed_jobs)
                + len(self.failed_jobs)
                + len(self.active_jobs)
                + preprocessing_queue_size
                + bed_conversion_queue_size
                + mgmt_queue_size
                + cnv_queue_size
                + target_queue_size
                + fusion_queue_size
                + classification_queue_size
                + slow_queue_size
            )
        else:
            total_expected = (
                len(self.completed_jobs)
                + len(self.failed_jobs)
                + len(self.active_jobs)
                + preprocessing_queue_size
                + bed_conversion_queue_size
                + analysis_queue_size
                + classification_queue_size
                + slow_queue_size
            )

        # If we have no jobs in progress but have completed jobs, use the completed count as the total
        if total_expected == 0 and len(self.completed_jobs) > 0:
            total_expected = len(self.completed_jobs)

        # Ensure total_expected is at least as large as the number of completed jobs
        if total_expected < len(self.completed_jobs):
            total_expected = len(self.completed_jobs)

        # Calculate total jobs that were actually enqueued (excluding skipped jobs)
        total_actual_jobs = self.total_jobs_enqueued

        # Prepare samples payload (copy to avoid race conditions)
        samples_payload: List[Dict[str, Any]] = []
        with self.progress_lock:
            for sid, info in self.samples_by_id.items():
                samples_payload.append(
                    {
                        "sample_id": info["sample_id"],
                        "active_jobs": info["active_jobs"],
                        "total_jobs": info["total_jobs"],
                        "completed_jobs": info["completed_jobs"],
                        "failed_jobs": info["failed_jobs"],
                        "job_types": list(info["job_types"]),
                        "last_seen": info["last_seen"],
                    }
                )

        return {
            "completed": len(self.completed_jobs),
            "failed": len(self.failed_jobs),
            "total_enqueued": self.total_jobs_enqueued,
            "total_skipped": self.total_jobs_skipped,
            "total_actual_jobs": total_actual_jobs,
            "total_processed": len(self.completed_jobs) + len(self.failed_jobs),
            "total_expected": total_expected,
            "active_jobs": len(self.active_jobs),
            "completed_by_type": self.completed_jobs_by_type.copy(),  # Make a copy to avoid threading issues
            "failed_by_type": self.failed_jobs_by_type.copy(),  # Make a copy to avoid threading issues
            "samples": samples_payload,
            "queue_sizes": {
                "preprocessing": preprocessing_queue_size,
                "bed_conversion": bed_conversion_queue_size,
                "mgmt": mgmt_queue_size,
                "cnv": cnv_queue_size,
                "target": target_queue_size,
                "fusion": fusion_queue_size,
                "analysis": analysis_queue_size,
                "classification": classification_queue_size,
                "slow": slow_queue_size,
            },
            "active_by_worker": active_by_worker,
            # Coalesce metrics (dispatch-time merging of small batches)
            "coalesce_events": int(getattr(self, "coalesce_events", 0) or 0),
            "coalesce_jobs_absorbed": int(
                getattr(self, "coalesce_jobs_absorbed", 0) or 0
            ),
            "coalesce_contexts_merged": int(
                getattr(self, "coalesce_contexts_merged", 0) or 0
            ),
        }

    # === Per-sample tracking helpers ===
    def _ensure_sample_entry(self, sample_id: str) -> Dict[str, Any]:
        with self.progress_lock:
            if sample_id not in self.samples_by_id:
                self.samples_by_id[sample_id] = {
                    "sample_id": sample_id,
                    "active_jobs": 0,
                    "total_jobs": 0,
                    "completed_jobs": 0,
                    "failed_jobs": 0,
                    "job_types": set(),
                    "last_seen": time.time(),
                }
            return self.samples_by_id[sample_id]

    def _on_sample_job_started(self, sample_id: str, job_type: str) -> None:
        entry = self._ensure_sample_entry(sample_id)
        with self.progress_lock:
            entry["active_jobs"] += 1
            entry["total_jobs"] += 1
            if job_type:
                entry["job_types"].add(job_type)
            entry["last_seen"] = time.time()

    def _on_sample_job_finished(
        self, sample_id: str, job_type: str, success: bool
    ) -> None:
        entry = self._ensure_sample_entry(sample_id)
        with self.progress_lock:
            if entry["active_jobs"] > 0:
                entry["active_jobs"] -= 1
            if success:
                entry["completed_jobs"] += 1
            else:
                entry["failed_jobs"] += 1
            if job_type:
                entry["job_types"].add(job_type)
            entry["last_seen"] = time.time()


# === File Watcher ===
class FileWatcher(FileSystemEventHandler):
    """File watcher that integrates with the workflow manager."""

    def __init__(
        self,
        watch_dir: str,
        preprocessor_func: Callable[[str], List[Job]],
        manager: WorkflowManager,
        target_panel: str,
        recursive: bool = True,
        patterns: Optional[List[str]] = None,
        ignore_patterns: Optional[List[str]] = None,
        verbose: bool = False,
        show_progress: bool = True,
    ):
        self.watch_dir = watch_dir
        self.preprocessor_func = preprocessor_func
        self.manager = manager
        self.recursive = recursive
        self.patterns = patterns or ["*"]
        self.ignore_patterns = ignore_patterns or []
        self.verbose = verbose
        self.show_progress = show_progress
        self.target_panel = target_panel
        self.observer = Observer()
        self.processed_files = set()

    def _should_process_file(self, filepath: str) -> bool:
        """Check if a file should be processed based on patterns."""
        path = Path(filepath)

        # Check if file matches any ignore patterns
        for pattern in self.ignore_patterns:
            if path.match(pattern):
                return False

        # Check if file matches any watch patterns
        for pattern in self.patterns:
            if path.match(pattern):
                return True

        return False

    def handle_file(self, filepath: str) -> None:
        """Process a new file through the workflow."""
        if filepath in self.processed_files:
            return

        self.processed_files.add(filepath)

        # Use logger for file watcher events
        import logging

        logger = logging.getLogger("robin.filewatcher")
        logger.info(f"Detected new file: {filepath}")

        try:
            # Check if preprocessor function accepts target_panel parameter
            import inspect
            sig = inspect.signature(self.preprocessor_func)
            if "target_panel" in sig.parameters:
                jobs = self.preprocessor_func(filepath, target_panel=self.target_panel)
            else:
                jobs = self.preprocessor_func(filepath)
            self.manager.enqueue_jobs(jobs)
            logger.debug(f"Queued {len(jobs)} jobs for {filepath}")
        except Exception as e:
            logger.error(f"Error processing {filepath}: {e}")

    def on_created(self, event) -> None:
        """Handle file creation events."""
        if not event.is_directory and self._should_process_file(event.src_path):
            self.handle_file(event.src_path)

    def on_modified(self, event) -> None:
        """Handle file modification events."""
        if not event.is_directory and self._should_process_file(event.src_path):
            self.handle_file(event.src_path)

    def on_moved(self, event) -> None:
        """Handle file move events (common with rsync default behavior)."""
        if not event.is_directory and self._should_process_file(event.dest_path):
            self.handle_file(event.dest_path)  # Use dest_path for the final location

    def start(self, process_existing: bool = True) -> None:
        """Start watching the directory."""
        if process_existing:
            self._process_existing_files()

        self.observer.schedule(self, self.watch_dir, recursive=self.recursive)
        self.observer.start()

        # Use logger for file watcher events
        import logging

        logger = logging.getLogger("robin.filewatcher")
        logger.info(f"Watching directory: {self.watch_dir}")
        logger.debug(f"Patterns: {self.patterns}")
        logger.debug(f"Ignore patterns: {self.ignore_patterns}")

    def _process_existing_files(self) -> None:
        """Process existing files that match the patterns."""
        import logging

        logger = logging.getLogger("robin.filewatcher")
        logger.info(f"Processing existing files in: {self.watch_dir}")

        existing_files = []

        # Find all existing files that match patterns
        for pattern in self.patterns:
            if pattern == "*":
                # Default pattern - all files
                if self.recursive:
                    pattern_files = list(Path(self.watch_dir).rglob("*"))
                else:
                    pattern_files = list(Path(self.watch_dir).glob("*"))
            else:
                # Simple pattern like *.bam
                if self.recursive:
                    pattern_files = list(Path(self.watch_dir).rglob(pattern))
                else:
                    pattern_files = list(Path(self.watch_dir).glob(pattern))

            # Filter to only files (not directories) and apply ignore patterns
            for file_path in pattern_files:
                if file_path.is_file():
                    # Check if file should be ignored
                    should_ignore = False
                    if self.ignore_patterns:
                        for ignore_pattern in self.ignore_patterns:
                            if file_path.match(ignore_pattern):
                                should_ignore = True
                                break

                    if not should_ignore:
                        existing_files.append(file_path)

        # Remove duplicates and sort
        existing_files = sorted(set(existing_files))

        if existing_files:
            logger.info(f"Found {len(existing_files)} existing file(s) to process:")
            for file_path in existing_files:
                logger.debug(f"  - {file_path}")

            # Process each existing file with progress bar
            with tqdm(
                total=len(existing_files),
                desc="Processing existing files",
                unit="file",
                disable=not self.show_progress,
            ) as pbar:
                for file_path in existing_files:
                    if self.show_progress:
                        pbar.set_postfix_str(f"Processing: {file_path.name}")
                    logger.debug(f"Processing existing file: {file_path}")
                    self.handle_file(str(file_path))
                    pbar.update(1)
        else:
            logger.info("No existing files found matching the patterns.")

    def stop(self, timeout: float = 30.0) -> bool:
        """Stop watching and shutdown the workflow manager."""
        if self.verbose:
            print("[Watcher] Stopping file watcher and workflow manager")

        # Stop the file observer
        self.observer.stop()
        self.observer.join()

        # Stop the workflow manager
        manager_stopped = self.manager.stop(timeout=timeout)

        if self.verbose:
            print("[Watcher] Stopped")

        return manager_stopped


# === Job Classifier ===
_job_id_counter = itertools.count(1000)


def default_file_classifier(filepath: str, workflow_plan: List[str], target_panel: str) -> List[Job]:
    """Default classifier that creates jobs for a file based on a workflow plan."""
    job_id = next(_job_id_counter)
    ctx = WorkflowContext(filepath)
    ctx.add_metadata("filename", os.path.basename(filepath))
    ctx.add_metadata("created", time.time())
    ctx.add_metadata("target_panel", target_panel)  # Add panel metadata

    # Try to get file size, but don't fail if file doesn't exist
    try:
        ctx.add_metadata("file_size", os.path.getsize(filepath))
    except (OSError, FileNotFoundError):
        ctx.add_metadata("file_size", 0)

    # Define job dependencies and triggers for parallel execution
    job_dependencies = {
        "bed_conversion": {"preprocessing"},
        "mgmt": {"preprocessing"},
        "cnv": {"preprocessing"},
        "target": {"preprocessing"},
        "fusion": {"preprocessing"},
        "sturgeon": {"bed_conversion"},
        "nanodx": {"bed_conversion"},
        "pannanodx": {"bed_conversion"},
        "random_forest": {"bed_conversion"},
    }

    # Create preprocessing job if it's the first step
    if workflow_plan and workflow_plan[0].endswith(":preprocessing"):
        # Determine what jobs this preprocessing can trigger
        triggers = set()
        for job_type, deps in job_dependencies.items():
            if "preprocessing" in deps:
                triggers.add(job_type)

        # Store the original workflow plan in context for reference
        ctx.add_metadata("original_workflow", workflow_plan)

        return [
            Job(
                job_id=job_id,
                job_type="preprocessing",
                context=ctx,
                origin="preprocessing",
                workflow=workflow_plan,
                step=0,
                dependencies=set(),
                triggers=triggers,
            )
        ]

    # Create first job in plan (skip preprocessing if not present)
    if len(workflow_plan) > 0:
        first_step = workflow_plan[0]
        if not isinstance(first_step, str):
            raise ValueError(
                f"Expected string in workflow step, got {type(first_step)}: {first_step}"
            )

        if ":" not in first_step:
            raise ValueError(
                f"Invalid workflow step format (missing ':'): {first_step}"
            )

        queue_type, job_type = first_step.split(":", 1)
        dependencies = job_dependencies.get(job_type, set())
        triggers = set()
        for trigger_job, trigger_deps in job_dependencies.items():
            if job_type in trigger_deps:
                triggers.add(trigger_job)

        # Store the original workflow plan in context for reference
        ctx.add_metadata("original_workflow", workflow_plan)

        return [
            Job(
                job_id=job_id,
                job_type=job_type,
                context=ctx,
                origin=queue_type,
                workflow=workflow_plan,
                step=0,
                dependencies=dependencies,
                triggers=triggers,
            )
        ]

    return []


# === Built-in Job Handlers ===
def echo_handler(job: Job) -> None:
    """Simple echo handler for testing."""
    time.sleep(0.1)
    ctx = job.context
    print(f"  [Echo] Processing {ctx.filepath}")
    ctx.add_metadata("echo_processed", True)
    ctx.add_result("echo", "echo_result")


def sleep_handler(job: Job) -> None:
    """Sleep handler for testing slow operations."""
    time.sleep(1.0)
    ctx = job.context
    print(f"  [Sleep] Processing {ctx.filepath}")
    ctx.add_metadata("sleep_processed", True)
    ctx.add_result("sleep", "sleep_result")


def command_handler(job: Job, command_template: str) -> None:
    """Handler that runs shell commands."""
    import subprocess

    ctx = job.context
    command = command_template.replace("{file}", ctx.filepath)

    try:
        result = subprocess.run(
            command,
            shell=True,
            capture_output=True,
            text=True,
            timeout=30,
        )

        ctx.add_result(
            job.job_type,
            {
                "command": command,
                "stdout": result.stdout,
                "stderr": result.stderr,
                "returncode": result.returncode,
            },
        )

        if result.returncode != 0:
            ctx.add_error(
                job.job_type, f"Command failed with return code {result.returncode}"
            )

    except subprocess.TimeoutExpired:
        ctx.add_error(job.job_type, "Command timed out")
    except Exception as e:
        ctx.add_error(job.job_type, f"Command error: {str(e)}")


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
            context.add_result(job_type, f"{job_type}_ok")
            
        except Exception as e:
            # Fail entire batch if any file fails
            error_msg = f"File {i+1}/{batch_size} failed: {str(e)}"
            for ctx in job.contexts:
                ctx.add_error(job_type, error_msg)
            raise  # Re-raise to fail the entire batch
    
    # Mark batch as completed
    for context in job.contexts:
        context.add_result(job_type, f"{job_type}_batch_completed")

def process_single_file_in_batch(context: WorkflowContext, job_type: str, index: int, total: int) -> None:
    """Process a single file within a batch - to be implemented per job type"""
    # This would call the existing single-file processing logic
    pass


# === Workflow Runner ===
class WorkflowRunner:
    """High-level interface for running workflows."""

    def __init__(
        self,
        target_panel: str,
        verbose: bool = False,
        analysis_workers: int = 1,
        use_separate_analysis_queues: bool = True,
        preprocessing_workers: int = 1,
        bed_workers: int = 1,
        reference: Optional[Path] = None,
        center: str = None,
        enable_batching: bool = True,
    ):
        self.manager = WorkflowManager(
            verbose=verbose,
            analysis_workers=analysis_workers,
            use_separate_analysis_queues=use_separate_analysis_queues,
            preprocessing_workers=preprocessing_workers,
            bed_conversion_workers=bed_workers,
            enable_batching=enable_batching,
        )
        self.verbose = verbose
        self.reference = reference
        self.center = center
        self.target_panel = target_panel

        # Log reference genome status

        # Log reference genome status
        if self.reference:
            print(f"[WorkflowRunner] Reference genome configured: {self.reference}")
        else:
            print(
                "[WorkflowRunner] No reference genome configured - SNP calling will not be available"
            )

        # Register default handlers
        self.manager.register_handler("preprocessing", "echo", echo_handler)
        self.manager.register_handler("slow", "sleep", sleep_handler)
        # Register IGV BAM handler on slow queue
        self.manager.register_handler("slow", "igv_bam", igv_bam_handler)
        # Register SNP analysis handler on slow queue
        self.manager.register_handler("slow", "snp_analysis", snp_analysis_handler)
        # Register target BAM finalization handler on slow queue
        self.manager.register_handler(
            "slow", "target_bam_finalize", target_bam_finalize_handler
        )

    def register_handler(
        self, queue_type: str, job_type: str, handler: Callable[[Job], None]
    ) -> None:
        """Register a custom job handler."""
        self.manager.register_handler(queue_type, job_type, handler)
    
    def register_batched_handler(
        self, job_type: str, handler: Callable[[BatchedJob], None]
    ) -> None:
        """Register a handler that can process batched jobs."""
        self.manager.register_batched_handler(job_type, handler)

    def register_command_handler(
        self, queue_type: str, job_type: str, command_template: str
    ) -> None:
        """Register a command handler that runs shell commands."""

        def handler(job: Job) -> None:
            command_handler(job, command_template)

        self.manager.register_handler(queue_type, job_type, handler)

    def submit_sample_job(
        self, sample_dir: str, job_type: str, sample_id: str = None
    ) -> bool:
        """
        Submit a job for an existing sample directory.

        This allows users to manually trigger specific job types for samples
        that have already been processed or need reprocessing.

        Args:
            sample_dir: Path to the sample directory
            job_type: Type of job to run (e.g., 'igv_bam')
            sample_id: Optional sample ID (defaults to directory name)

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
                    "bam_metadata": {"sample_id": sample_id},
                },
            )

            # Create a job
            job = Job(
                job_id=next(_job_id_counter),
                job_type=job_type,
                context=context,
                origin="slow",  # IGV BAM jobs go to slow queue
                workflow=[f"slow:{job_type}"],
            )

            # Submit the job to the appropriate queue
            self.manager.enqueue_jobs([job])

            if self.verbose:
                print(
                    f"[WorkflowRunner] Submitted {job_type} job for sample {sample_id}"
                )

            return True

        except Exception as e:
            if self.verbose:
                print(
                    f"[WorkflowRunner] Failed to submit {job_type} job for sample {sample_id}: {e}"
                )
            return False

    def submit_snp_analysis_job(
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
        try:
            if sample_id is None:
                sample_id = Path(sample_dir).name

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
                "bam_metadata": {"sample_id": sample_id},
                "threads": threads,
                "force_regenerate": force_regenerate,
            }

            if target_panel:
                metadata["target_panel"] = target_panel

            # Add reference genome if provided
            if reference:
                metadata["reference"] = reference

            context = WorkflowContext(filepath=sample_dir, metadata=metadata)

            # Create a job
            job = Job(
                job_id=next(_job_id_counter),
                job_type="snp_analysis",
                context=context,
                origin="slow",  # SNP analysis jobs go to slow queue
                workflow=["slow:snp_analysis"],
            )

            # Submit the job to the appropriate queue
            self.manager.enqueue_jobs([job])

            if self.verbose:
                print(
                    f"[WorkflowRunner] Submitted SNP analysis job for sample {sample_id}"
                )
                if reference:
                    print(f"[WorkflowRunner] Using reference genome: {reference}")

            return True

        except Exception as e:
            if self.verbose:
                print(
                    f"[WorkflowRunner] Failed to submit SNP analysis job for sample {sample_id}: {e}"
                )
            return False

    def submit_target_bam_finalize_job(
        self,
        sample_dir: str,
        sample_id: str = None,
        target_panel: str = None,
    ) -> bool:
        """
        Submit a target BAM finalization job for an existing sample directory.
        """
        try:
            if sample_id is None:
                sample_id = Path(sample_dir).name

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

            context = WorkflowContext(filepath=sample_dir, metadata=metadata)

            job = Job(
                job_id=next(_job_id_counter),
                job_type="target_bam_finalize",
                context=context,
                origin="slow",
                workflow=["slow:target_bam_finalize"],
            )

            self.manager.enqueue_jobs([job])

            if self.verbose:
                print(
                    f"[WorkflowRunner] Submitted target BAM finalization job for sample {sample_id}"
                )

            return True
        except Exception as e:
            if self.verbose:
                print(
                    f"[WorkflowRunner] Failed to submit target BAM finalization job for sample {sample_id}: {e}"
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

    def run_workflow(
        self,
        watch_dir: str,
        workflow_plan: List[str],
        recursive: bool = True,
        patterns: Optional[List[str]] = None,
        ignore_patterns: Optional[List[str]] = None,
        classifier_func: Optional[Callable] = None,
        process_existing: bool = True,
        show_progress: bool = True,
    ) -> None:
        """Run a complete workflow."""

        if classifier_func is None:
            classifier_func = lambda filepath: default_file_classifier(
                filepath, workflow_plan, self.target_panel
            )

        watcher = FileWatcher(
            watch_dir=watch_dir,
            preprocessor_func=classifier_func,
            manager=self.manager,
            target_panel=self.target_panel,
            recursive=recursive,
            patterns=patterns,
            ignore_patterns=ignore_patterns,
            verbose=self.verbose,
            show_progress=show_progress,
        )

        # Start worker threads first to ensure queues are consumed immediately
        self.manager.run()

        # Now start watching files (which enqueues jobs)
        watcher.start(process_existing=process_existing)

        # Start progress monitoring if enabled
        progress_thread = None
        if show_progress:
            progress_thread = threading.Thread(
                target=self._monitor_progress, args=(watcher,), daemon=True
            )
            progress_thread.start()

        try:
            # Keep the main thread alive while workers and watcher run
            while self.manager.is_running():
                time.sleep(0.5)
        except KeyboardInterrupt:
            print("\n[SHUTDOWN] Interrupted by user (Ctrl-C)")
            print("[SHUTDOWN] Initiating graceful shutdown...")
            if self.verbose:
                print("[WorkflowRunner] Shutdown requested.")
            print("[WorkflowRunner] Shutdown requested.")
            # Stop the watcher and workflow manager
            print("[SHUTDOWN] Stopping file watcher...")
            graceful_shutdown = watcher.stop(timeout=1.0)
            if graceful_shutdown:
                print("[SHUTDOWN] File watcher stopped gracefully")
            else:
                print("[SHUTDOWN] Warning: File watcher may not have stopped gracefully")

            if self.verbose:
                if graceful_shutdown:
                    print("[WorkflowRunner] Graceful shutdown completed.")
                else:
                    print(
                        "[WorkflowRunner] Warning: Some workers may not have stopped gracefully."
                    )
                stats = self.manager.get_stats()
                print(f"[WorkflowRunner] Final stats: {stats}")
                print(
                    f"[WorkflowRunner] Final completion by type: {stats.get('completed_by_type', {})}"
                )
                print(
                    f"[WorkflowRunner] Final failed by type: {stats.get('failed_by_type', {})}"
                )

    def _monitor_progress(self, watcher) -> None:
        """Monitor and display worker progress in real-time."""
        import time
        if _should_use_rich_progress():
            self._monitor_progress_rich()
            return

        # Create progress bars for each worker type (tqdm fallback)
        preprocessing_pbar = tqdm(
            desc="Preprocessing", unit="jobs", position=0, leave=True
        )
        bed_conversion_pbar = tqdm(
            desc="Bed Conversion", unit="jobs", position=1, leave=True
        )
        mgmt_pbar = tqdm(desc="MGMT", unit="jobs", position=2, leave=True)
        cnv_pbar = tqdm(desc="CNV", unit="jobs", position=3, leave=True)
        target_pbar = tqdm(desc="Target", unit="jobs", position=4, leave=True)
        fusion_pbar = tqdm(desc="Fusion", unit="jobs", position=5, leave=True)
        classification_pbar = tqdm(
            desc="Classification", unit="jobs", position=6, leave=True
        )
        slow_pbar = tqdm(desc="Slow", unit="jobs", position=7, leave=True)

        # Create overall progress bar
        overall_pbar = tqdm(
            desc="Overall Progress", unit="jobs", position=8, leave=True
        )

        progress_bars = {
            "preprocessing": preprocessing_pbar,
            "bed_conversion": bed_conversion_pbar,
            "mgmt": mgmt_pbar,
            "cnv": cnv_pbar,
            "target": target_pbar,
            "fusion": fusion_pbar,
            "classification": classification_pbar,
            "slow": slow_pbar,
        }

        try:
            while self.manager.is_running():
                stats = self.manager.get_stats()

                # Update overall progress using total actual jobs (excluding skipped)
                total_processed = stats["total_processed"]
                total_actual_jobs = stats["total_actual_jobs"]

                if total_actual_jobs > 0:
                    overall_pbar.total = total_actual_jobs
                    overall_pbar.n = total_processed
                    overall_pbar.set_postfix_str(
                        f"Active: {stats['active_jobs']} | "
                        f"Completed: {stats['completed']} | "
                        f"Failed: {stats['failed']} | "
                        f"Skipped: {stats['total_skipped']} | "
                        f"Total: {total_actual_jobs}"
                    )

                # Update queue-specific progress bars
                for queue_name, pbar in progress_bars.items():
                    queue_size = stats["queue_sizes"][queue_name]

                    # Count active jobs for this queue type
                    active_in_queue = 0
                    active_jobs_info = []

                    for worker_name, jobs in stats["active_by_worker"].items():
                        for job in jobs:
                            # Check if this worker belongs to the current queue
                            if (
                                queue_name == "preprocessing"
                                and worker_name.startswith("PreprocessingWorker")
                            ):
                                active_in_queue += 1
                                filename = job["filepath"].split("/")[-1]
                                duration = int(job["duration"])
                                active_jobs_info.append(
                                    f"{job['job_type']}:{filename}({duration}s)"
                                )
                            elif (
                                queue_name == "bed_conversion"
                                and worker_name.startswith("BedConversionWorker")
                            ):
                                active_in_queue += 1
                                filename = job["filepath"].split("/")[-1]
                                duration = int(job["duration"])
                                active_jobs_info.append(
                                    f"{job['job_type']}:{filename}({duration}s)"
                                )
                            elif queue_name == "mgmt" and worker_name.startswith(
                                "MGMTWorker"
                            ):
                                active_in_queue += 1
                                filename = job["filepath"].split("/")[-1]
                                duration = int(job["duration"])
                                active_jobs_info.append(
                                    f"{job['job_type']}:{filename}({duration}s)"
                                )
                            elif queue_name == "cnv" and worker_name.startswith(
                                "CNVWorker"
                            ):
                                active_in_queue += 1
                                filename = job["filepath"].split("/")[-1]
                                duration = int(job["duration"])
                                active_jobs_info.append(
                                    f"{job['job_type']}:{filename}({duration}s)"
                                )
                            elif queue_name == "target" and worker_name.startswith(
                                "TargetWorker"
                            ):
                                active_in_queue += 1
                                filename = job["filepath"].split("/")[-1]
                                duration = int(job["duration"])
                                active_jobs_info.append(
                                    f"{job['job_type']}:{filename}({duration}s)"
                                )
                            elif queue_name == "fusion" and worker_name.startswith(
                                "FusionWorker"
                            ):
                                active_in_queue += 1
                                filename = job["filepath"].split("/")[-1]
                                duration = int(job["duration"])
                                active_jobs_info.append(
                                    f"{job['job_type']}:{filename}({duration}s)"
                                )
                            elif (
                                queue_name == "classification"
                                and worker_name.startswith("ClassificationWorker")
                            ):
                                active_in_queue += 1
                                filename = job["filepath"].split("/")[-1]
                                duration = int(job["duration"])
                                active_jobs_info.append(
                                    f"{job['job_type']}:{filename}({duration}s)"
                                )
                            elif queue_name == "slow" and worker_name.startswith(
                                "SlowWorker"
                            ):
                                active_in_queue += 1
                                filename = job["filepath"].split("/")[-1]
                                duration = int(job["duration"])
                                active_jobs_info.append(
                                    f"{job['job_type']}:{filename}({duration}s)"
                                )

                    # Calculate completed jobs for this queue type (both successful and failed)
                    completed_in_queue = 0
                    for job_type, count in stats["completed_by_type"].items():
                        # Get queue type from the manager's job queue mapping
                        queue_type = self.manager.job_queue_mapping.get(
                            job_type, "slow"
                        )
                        if queue_type == queue_name:
                            completed_in_queue += count

                    # Add failed jobs to the completion count
                    for job_type, count in stats.get("failed_by_type", {}).items():
                        # Get queue type from the manager's job queue mapping
                        queue_type = self.manager.job_queue_mapping.get(
                            job_type, "slow"
                        )
                        if queue_type == queue_name:
                            completed_in_queue += count

                    # Calculate total jobs for this queue (completed + active + queued)
                    total_for_queue = completed_in_queue + active_in_queue + queue_size

                    # Update the progress bar with total and current values
                    if total_for_queue > 0:
                        pbar.total = total_for_queue
                        pbar.n = completed_in_queue

                    # Update the progress bar description with queue info
                    pbar.set_description(
                        f"{queue_name.title()} (Q:{queue_size} A:{active_in_queue} C:{completed_in_queue})"
                    )

                    # Update active jobs info
                    if active_jobs_info:
                        pbar.set_postfix_str(
                            " | ".join(active_jobs_info[:2])
                        )  # Show first 2 active jobs

                # Check if we should stop monitoring
                if not self.manager.is_running():
                    break

                time.sleep(1.0)  # Update every second

        except KeyboardInterrupt:
            pass
        finally:
            # Close all progress bars
            for pbar in progress_bars.values():
                pbar.close()
            overall_pbar.close()

    def _monitor_progress_rich(self) -> None:
        """Monitor and display worker progress using rich progress bars."""
        import time

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

        queue_order = [
            "preprocessing",
            "bed_conversion",
            "mgmt",
            "cnv",
            "target",
            "fusion",
            "classification",
            "slow",
        ]
        task_ids = {
            queue: progress.add_task(queue.title(), total=None, detail="")
            for queue in queue_order
        }
        overall_task_id = progress.add_task("Overall Progress", total=None, detail="")

        progress.start()
        try:
            while self.manager.is_running():
                stats = self.manager.get_stats()

                total_processed = stats["total_processed"]
                total_actual_jobs = stats["total_actual_jobs"]
                overall_total = total_actual_jobs if total_actual_jobs > 0 else None

                progress.update(
                    overall_task_id,
                    total=overall_total,
                    completed=total_processed,
                    detail=(
                        f"Done:{total_processed}/{total_actual_jobs} | "
                        f"Act:{stats['active_jobs']} | "
                        f"Fail:{stats['failed']} | "
                        f"Skip:{stats['total_skipped']}"
                    ),
                )

                for queue_name in queue_order:
                    queue_size = stats["queue_sizes"][queue_name]
                    active_in_queue = 0
                    active_jobs_info = []

                    for worker_name, jobs in stats["active_by_worker"].items():
                        for job in jobs:
                            filename = job["filepath"].split("/")[-1]
                            duration = int(job["duration"])
                            matches_queue = False
                            if (
                                queue_name == "preprocessing"
                                and worker_name.startswith("PreprocessingWorker")
                            ):
                                active_in_queue += 1
                                matches_queue = True
                            elif (
                                queue_name == "bed_conversion"
                                and worker_name.startswith("BedConversionWorker")
                            ):
                                active_in_queue += 1
                                matches_queue = True
                            elif queue_name == "mgmt" and worker_name.startswith(
                                "MGMTWorker"
                            ):
                                active_in_queue += 1
                                matches_queue = True
                            elif queue_name == "cnv" and worker_name.startswith(
                                "CNVWorker"
                            ):
                                active_in_queue += 1
                                matches_queue = True
                            elif queue_name == "target" and worker_name.startswith(
                                "TargetWorker"
                            ):
                                active_in_queue += 1
                                matches_queue = True
                            elif queue_name == "fusion" and worker_name.startswith(
                                "FusionWorker"
                            ):
                                active_in_queue += 1
                                matches_queue = True
                            elif (
                                queue_name == "classification"
                                and worker_name.startswith("ClassificationWorker")
                            ):
                                active_in_queue += 1
                                matches_queue = True
                            elif queue_name == "slow" and worker_name.startswith(
                                "SlowWorker"
                            ):
                                active_in_queue += 1
                                matches_queue = True

                            if matches_queue and len(active_jobs_info) < 2:
                                active_jobs_info.append(
                                    f"{job['job_type']}:{filename}({duration}s)"
                                )

                    completed_in_queue = 0
                    for job_type, count in stats["completed_by_type"].items():
                        queue_type = self.manager.job_queue_mapping.get(
                            job_type, "slow"
                        )
                        if queue_type == queue_name:
                            completed_in_queue += count
                    for job_type, count in stats.get("failed_by_type", {}).items():
                        queue_type = self.manager.job_queue_mapping.get(
                            job_type, "slow"
                        )
                        if queue_type == queue_name:
                            completed_in_queue += count

                    total_for_queue = completed_in_queue + active_in_queue + queue_size
                    total_for_queue_display = (
                        str(total_for_queue) if total_for_queue > 0 else "-"
                    )
                    completed_display = str(completed_in_queue)
                    detail_parts = [f"{completed_display}/{total_for_queue_display}"]
                    if active_jobs_info:
                        detail_parts.append(" | ".join(active_jobs_info))
                    progress.update(
                        task_ids[queue_name],
                        total=total_for_queue if total_for_queue > 0 else None,
                        completed=completed_in_queue,
                        detail=" | ".join(detail_parts),
                    )

                if not self.manager.is_running():
                    break

                time.sleep(1.0)
        finally:
            progress.stop()


def _should_use_rich_progress() -> bool:
    if not _RICH_AVAILABLE:
        return False
    preference = os.environ.get("ROBIN_PROGRESS", "").strip().lower()
    if preference in {"tqdm", "plain", "off", "0", "false", "no"}:
        return False
    if preference in {"rich", "on", "1", "true", "yes"}:
        return True
    return True
