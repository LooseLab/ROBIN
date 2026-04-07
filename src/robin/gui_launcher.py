"""
GUI launcher for robin workflow monitoring.

This module provides a clean interface for launching the workflow GUI
that runs in a separate thread but is completely isolated from
workflow execution to avoid any blocking.
"""

# Suppress pkg_resources deprecation warnings from sorted_nearest
import warnings
warnings.filterwarnings(
    "ignore", message="pkg_resources is deprecated", category=UserWarning
)
# Suppress matplotlib tight_layout warnings
warnings.filterwarnings(
    "ignore", message="The figure layout has changed to tight", category=UserWarning
)

import asyncio
import hashlib
import logging
import queue
import threading
import time
from collections import deque
import csv
from datetime import datetime
from robin.analysis.master_csv_manager import MasterCSVManager

from typing import Optional, Dict, Any, List, Set
from pathlib import Path
from dataclasses import dataclass, field, asdict
from enum import Enum
import os
import secrets
import sys
import tempfile
import zipfile
import json
import pickle
import getpass
from urllib.parse import quote

from robin.gui import theme, images

from robin.gui.components.news_feed import NewsFeed

from robin.reporting.report import create_pdf
from robin.reporting.sections.disclaimer_text import EXTENDED_DISCLAIMER_TEXT


# Files that indicate an analysis step is complete.
COMPLETION_JOB_PATTERNS: Dict[str, List[str]] = {
    "fusion": [
        "fusion_candidates_master_processed.pkl",
        "fusion_candidates_all_processed.pkl",
        "fusion_summary.csv",
        "fusion_results.csv",
        "sv_count.txt",
    ],
    "cnv": [
        "CNV_dict.npy",
        "XYestimate.pkl",
        "cnv_results.csv",
        "cnv_summary.txt",
        "cnv_analysis_counter.txt",
        "cnv_analysis_results.pkl",
    ],
    "mgmt": ["final_mgmt.csv", "*_mgmt.csv"],
    "target": ["coverage_main.csv", "bed_coverage_main.csv"],
    "sturgeon": ["sturgeon_scores.csv", "sturgeon_results.csv", "sturgeon_summary.csv"],
    "nanodx": ["NanoDX_scores.csv", "nanodx_results.csv", "nanodx_summary.csv"],
    "pannanodx": ["PanNanoDX_scores.csv", "pannanodx_results.csv", "pannanodx_summary.csv"],
    "random_forest": [
        "random_forest_scores.csv",
        "random_forest_results.csv",
        "random_forest_summary.csv",
    ],
}

# Manifest written when a sample ID is generated from Test ID / name / DOB (for matching sequencer data).
SAMPLE_IDENTIFIER_MANIFEST_FILENAME = "sample_identifier_manifest.json"


def _get_test_id_from_manifest(sample_dir: Path) -> str:
    """Read test_id from sample_identifier_manifest.json in the sample directory if present."""
    manifest_path = sample_dir / SAMPLE_IDENTIFIER_MANIFEST_FILENAME
    if not manifest_path.exists():
        return ""
    try:
        with open(manifest_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        return str(data.get("test_id", "") or "").strip()
    except Exception:
        return ""

# Salt for deriving encryption key from DOB (fixed so the same DOB always produces the same key).
_IDENTIFIER_MANIFEST_KEY_SALT = b"robin_sample_manifest_v1"


def _derive_key_from_dob(dob: str) -> bytes:
    """Derive a Fernet key from the date of birth for encrypting manifest fields."""
    import base64
    from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC
    from cryptography.hazmat.primitives import hashes

    kdf = PBKDF2HMAC(
        algorithm=hashes.SHA256(),
        length=32,
        salt=_IDENTIFIER_MANIFEST_KEY_SALT,
        iterations=100_000,
    )
    key_bytes = kdf.derive(dob.encode("utf-8"))
    return base64.urlsafe_b64encode(key_bytes)


def _encrypt_identifier_manifest_field(plaintext: str, dob: str) -> str:
    """Encrypt a string using a key derived from DOB; returns base64-encoded ciphertext."""
    from cryptography.fernet import Fernet

    key = _derive_key_from_dob(dob)
    f = Fernet(key)
    ciphertext = f.encrypt(plaintext.encode("utf-8"))
    return ciphertext.decode("ascii")


def _decrypt_identifier_manifest_field(ciphertext_b64: str, dob: str) -> str:
    """Decrypt a base64-encoded ciphertext using a key derived from DOB."""
    from cryptography.fernet import Fernet

    key = _derive_key_from_dob(dob)
    f = Fernet(key)
    plaintext = f.decrypt(ciphertext_b64.encode("ascii"))
    return plaintext.decode("utf-8")


def _load_manifest_encrypted_fields(sample_dir: Optional[Path]) -> Optional[Dict[str, str]]:
    """Load encrypted first_name, last_name, dob, nhs_number (hospital number) from sample_identifier_manifest.json if present."""
    if not sample_dir or not sample_dir.exists():
        return None
    manifest_path = sample_dir / SAMPLE_IDENTIFIER_MANIFEST_FILENAME
    if not manifest_path.exists():
        return None
    try:
        with open(manifest_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        fn = data.get("first_name")
        ln = data.get("last_name")
        dob_enc = data.get("dob")
        if fn and ln and dob_enc:
            out: Dict[str, str] = {"first_name": fn, "last_name": ln, "dob": dob_enc}
            nhs_enc = data.get("nhs_number")
            if nhs_enc:
                out["nhs_number"] = nhs_enc
            return out
    except Exception:
        pass
    return None


try:
    from nicegui import ui, app
except ImportError:
    ui = None
    app = None

try:
    from fastapi import Request
    from fastapi.responses import RedirectResponse
    from starlette.middleware.base import BaseHTTPMiddleware
except ImportError:  # pragma: no cover
    Request = None
    RedirectResponse = None
    BaseHTTPMiddleware = object

try:
    from argon2 import PasswordHasher
    from argon2.exceptions import VerifyMismatchError, InvalidHashError
except ImportError:  # pragma: no cover
    PasswordHasher = None
    VerifyMismatchError = Exception
    InvalidHashError = Exception


def _get_gui_password_hash_path() -> Path:
    """Return path to the file where the GUI password hash is stored."""
    if os.name == "nt":
        base = Path(os.environ.get("APPDATA", os.path.expanduser("~")))
    else:
        base = Path(os.environ.get("XDG_CONFIG_HOME", os.path.expanduser("~/.config")))
    return base / "robin" / "gui_password_hash"


def ensure_gui_password_set() -> bool:
    """Prompt for GUI password at startup: set (twice) if no hash exists, else verify once.

    Uses getpass so the password is not echoed. Returns True if the GUI may start,
    False otherwise (caller should exit). Requires a TTY to set a new password.
    """
    if PasswordHasher is None:
        logging.error("argon2-cffi is required for GUI password hashing. Install it with: pip install argon2-cffi")
        return False

    path = _get_gui_password_hash_path()
    hasher = PasswordHasher()

    if path.exists():
        try:
            stored = path.read_text().strip()
        except OSError as e:
            logging.error("Could not read GUI password hash file: %s", e)
            return False
        if not stored:
            logging.error("GUI password hash file is empty. Delete it and run again to set a new password.")
            return False
        if sys.stdin.isatty():
            try:
                pwd = getpass.getpass("GUI password: ")
            except (EOFError, KeyboardInterrupt):
                return False
            try:
                hasher.verify(stored, pwd)
            except VerifyMismatchError:
                print("Invalid password.", file=sys.stderr)
                return False
        return True

    if not sys.stdin.isatty():
        print(
            "No GUI password has been set. Run ROBIN from a terminal to set the password (you will be prompted twice).",
            file=sys.stderr,
        )
        return False

    try:
        pwd1 = getpass.getpass("Set GUI password: ")
        pwd2 = getpass.getpass("Confirm GUI password: ")
    except (EOFError, KeyboardInterrupt):
        return False

    if pwd1 != pwd2:
        print("Passwords do not match.", file=sys.stderr)
        return False
    if not pwd1:
        print("Password cannot be empty.", file=sys.stderr)
        return False

    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(hasher.hash(pwd1), encoding="utf-8")
    except OSError as e:
        logging.error("Could not write GUI password hash file: %s", e)
        return False
    return True


def set_gui_password_interactive() -> bool:
    """Interactive prompt to set or replace the GUI password (e.g. from `robin password set`).

    Uses getpass for all password input; nothing is echoed. If a password already exists,
    asks whether to replace it (no existing password required). Returns True on success.
    """
    if PasswordHasher is None:
        print(
            "argon2-cffi is required. Install it with: pip install argon2-cffi",
            file=sys.stderr,
        )
        return False

    console = None
    rich_confirm = None
    try:  # Use rich if available for nicer prompts/output
        from rich.console import Console  # type: ignore
        from rich.prompt import Confirm as RichConfirm  # type: ignore

        console = Console()
        rich_confirm = RichConfirm
    except Exception:  # pragma: no cover
        console = None
        rich_confirm = None

    def _echo(message: str, style: str = "") -> None:
        if console is not None and style:
            console.print(f"[{style}]{message}[/{style}]")
        elif console is not None:
            console.print(message)
        else:
            print(message)

    path = _get_gui_password_hash_path()
    hasher = PasswordHasher()

    if path.exists():
        try:
            if rich_confirm is not None:
                replace = rich_confirm.ask(
                    "A password is already set. Replace it?", default=False
                )
            else:
                reply = input("A password is already set. Replace it? [y/N]: ").strip().lower()
                replace = reply in ("y", "yes")
        except (EOFError, KeyboardInterrupt):
            return False
        if not replace:
            return False

    try:
        pwd1 = getpass.getpass("New GUI password: ")
        pwd2 = getpass.getpass("Confirm new GUI password: ")
    except (EOFError, KeyboardInterrupt):
        return False

    if pwd1 != pwd2:
        _echo("Passwords do not match.", style="bold red")
        return False
    if not pwd1:
        _echo("Password cannot be empty.", style="bold red")
        return False

    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(hasher.hash(pwd1), encoding="utf-8")
    except OSError as e:
        _echo(f"Could not write password file: {e}", style="bold red")
        return False

    _echo("GUI password updated successfully.", style="bold green")
    return True


class UpdateType(Enum):
    """Types of updates that can be sent to the GUI."""

    WORKFLOW_STATUS = "workflow_status"
    JOB_UPDATE = "job_update"
    QUEUE_UPDATE = "queue_update"
    LOG_MESSAGE = "log_message"
    PROGRESS_UPDATE = "progress_update"
    ERROR_UPDATE = "error_update"
    SAMPLES_UPDATE = "samples_update"
    WARNING_NOTIFICATION = "warning_notification"


@dataclass
class GUIUpdate:
    """A single update message for the GUI."""

    update_type: UpdateType
    timestamp: float
    data: Dict[str, Any]
    priority: int = 0  # Higher priority updates are processed first


@dataclass
class SampleRecord:
    """Master record for a single sample with all tracked information."""

    sample_id: str
    test_id: str = ""  # From sample_identifier_manifest.json if present; else empty
    origin: str = "Live"  # Live, Pre-existing, or Complete
    run_start: str = ""
    device: str = ""
    flowcell: str = ""
    active_jobs: int = 0
    pending_jobs: int = 0
    total_jobs: int = 0
    completed_jobs: int = 0
    failed_jobs: int = 0
    job_types: str = ""
    last_seen: str = ""  # Formatted timestamp string
    _last_seen_raw: float = 0.0  # Raw timestamp for sorting
    _dirty: bool = False  # Flag to indicate record needs UI update
    _file_mtime: float = 0.0  # master.csv modification time
    files_seen: int = 0  # Total files seen/processed
    files_processed: int = 0  # Files completely processed through all analysis steps
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for table rows."""
        return {
            "sample_id": self.sample_id,
            "test_id": self.test_id,
            "origin": self.origin,
            "run_start": self.run_start,
            "device": self.device,
            "flowcell": self.flowcell,
            "active_jobs": self.active_jobs,
            "pending_jobs": self.pending_jobs,
            "total_jobs": self.total_jobs,
            "completed_jobs": self.completed_jobs,
            "failed_jobs": self.failed_jobs,
            "job_types": self.job_types,
            "last_seen": self.last_seen,
            "files_seen": self.files_seen,
            "files_processed": self.files_processed,
            "file_progress": self.files_processed / self.files_seen if self.files_seen > 0 else 0.0,
            "actions": "View",
            "_last_seen_raw": self._last_seen_raw,
        }


class GUILauncher:
    """Launcher for the robin workflow GUI using isolated threading with message queue."""

    def __init__(self, host: str = "0.0.0.0", port: int = 8081, reload: bool = False):
        self.host = host
        self.port = port
        self.gui_thread = None
        self.is_running = False
        self.workflow_runner = None
        self.workflow_steps = []
        self.monitored_directory = ""
        self.reload = reload
        self.center = None  # Center ID for the analysis
        self._auth_middleware_registered = False
        # /robin_dark_mode: session cookie for theme (must work on /login before auth).
        self._unrestricted_page_routes = {"/login", "/robin_dark_mode"}
        self._password_hash: Optional[str] = None  # cached after first read

        # Sample status transition timeout (in seconds)
        # Change to 60 for testing (1 minute), 3600 for production (60 minutes)
        self.completion_timeout_seconds = 60  # Set to 3600 for production
        
        # Message queue for non-blocking communication
        self.update_queue = queue.PriorityQueue()
        self.gui_ready = threading.Event()
        self.shutdown_event = threading.Event()

        self.news_feed = None
        # GUI update thread
        self.update_thread = (
            None  # not used anymore; updates processed on UI thread via timer
        )

        # Debug counters
        self.total_updates_enqueued = 0
        self.total_updates_processed = 0
        # Monotonic tiebreaker for priority queue to avoid tuple comparison of GUIUpdate
        self._update_seq: int = 0
        
        # Adaptive throttling and update coalescing
        self._last_update_times: Dict[UpdateType, float] = {}
        self._last_queue_size_logged: int = 0
        self._last_log_time: float = 0.0

        # Runtime state
        self._start_time: Optional[float] = None
        self._is_running: bool = False
        
        # Track which samples have had target.bam finalization triggered
        self._finalized_samples: set = set()

        # SNP analysis concurrency control
        self._snp_analysis_lock = threading.Lock()
        self._snp_analysis_running_sample: Optional[str] = None

        # Log ring buffer (last 1000 entries)
        self._log_buffer = deque(maxlen=1000)

        # Cached samples data for persistence across navigation
        self._last_samples_rows: List[Dict[str, Any]] = []
        self._current_sample_id: Optional[str] = None
        self._selected_sample_id: Optional[str] = None
        self._known_sample_ids: set[str] = set()
        self._preexisting_sample_ids: set[str] = set()
        self._preexisting_scanned: bool = False
        self._pending_samples_data: Optional[Dict[str, Any]] = None
        
        # Caching for samples data
        self._last_cache_time: float = 0.0
        self._cache_duration: float = 30.0  # Cache for 30 seconds
        
        # Master record system - centralized source of truth for samples
        self._samples_master_record: Dict[str, SampleRecord] = {}
        # Per-GUI-process snapshot of on-disk totals from master.csv at first use per sample.
        # Coordinator stats are session-only; merged display = baseline + session counters.
        self._samples_job_baseline: Dict[str, Dict[str, int]] = {}
        # RLock: _merge_job_types_with_persisted and other helpers acquire this while callers
        # (e.g. _update_master_record_from_workflow) may already hold it — plain Lock deadlocks.
        self._samples_record_lock = threading.RLock()
        self._master_record_cache_file: Optional[Path] = None  # Set when monitored_directory is known
        self._background_scan_interval: float = 10.0  # Scan every 10 seconds
        self._background_scan_in_progress: bool = False
        # Sequential bulk SNP (one sample at a time; guard concurrent runs)
        self._bulk_snp_worker_running: bool = False
        # Per-sample manual pipeline status for finalize/SNP visibility in table.
        self._sample_pipeline_status: Dict[str, Dict[str, Any]] = {}

        # Sequential bulk MNP-Flex (one sample at a time; guard concurrent runs)
        self._bulk_mnpflex_worker_running: bool = False

        # MGMT per-sample cache (latest count seen)
        self._mgmt_state: Dict[str, Dict[str, Any]] = {}
        # Coverage per-sample cache (file mtimes, computed metrics)
        self._coverage_state: Dict[str, Dict[str, Any]] = {}
        
        # Progress notification event
        from nicegui import Event
        self.progress_notification_event = Event[Dict[str, Any]]()
        # CNV per-sample cache
        self._cnv_state: Dict[str, Dict[str, Any]] = {}
        # Cache last seen queue status so we can populate immediately on page creation
        self._last_queue_status: Dict[str, Any] = {}

    def _is_authenticated(self) -> bool:
        """Return True when the current browser session is authenticated for this server run.

        Requires both authenticated flag and matching _auth_generation so that persisted
        user storage from a previous server run (or before a password reset) is invalid.
        """
        if not app:
            return False
        gen = app.storage.general.get("_auth_generation")
        return bool(
            app.storage.user.get("authenticated", False)
            and gen is not None
            and app.storage.user.get("_auth_generation") == gen
        )

    def _register_auth_middleware(self) -> None:
        """Register request middleware that protects all GUI pages."""
        if self._auth_middleware_registered:
            return
        if app is None or BaseHTTPMiddleware is object:
            return

        # Invalidate any persisted auth from previous runs: each server run gets a new token.
        try:
            app.storage.general["_auth_generation"] = secrets.token_hex(16)
        except Exception:
            pass

        unrestricted_page_routes = set(self._unrestricted_page_routes)

        @app.add_middleware
        class AuthMiddleware(BaseHTTPMiddleware):
            async def dispatch(self, request: Request, call_next):
                # Store request host for theme (e.g. to show Quit only on localhost)
                try:
                    url = request.url
                    host = getattr(url, "hostname", None)
                    if not host and getattr(url, "host", None):
                        host = str(url.host).split(":")[0]
                    app.storage.user["_request_host"] = (host or "").strip()
                except Exception:
                    try:
                        app.storage.user["_request_host"] = ""
                    except Exception:
                        pass
                path = request.url.path
                if (
                    path.startswith("/_nicegui")
                    or path in unrestricted_page_routes
                    or path.startswith("/api/v1/")
                ):
                    return await call_next(request)
                gen = app.storage.general.get("_auth_generation")
                if not app.storage.user.get("authenticated", False) or gen is None or app.storage.user.get("_auth_generation") != gen:
                    requested_path = path
                    if request.url.query:
                        requested_path = f"{requested_path}?{request.url.query}"
                    redirect_to = quote(requested_path, safe="/?=&")
                    return RedirectResponse(url=f"/login?redirect_to={redirect_to}")
                return await call_next(request)

        self._auth_middleware_registered = True

    def _get_password_hash(self) -> Optional[str]:
        """Read and cache the stored password hash from the config file."""
        if self._password_hash is not None:
            return self._password_hash
        path = _get_gui_password_hash_path()
        if not path.exists():
            return None
        try:
            self._password_hash = path.read_text(encoding="utf-8").strip()
            return self._password_hash if self._password_hash else None
        except OSError:
            return None

    def _verify_password(self, candidate: str) -> bool:
        """Verify the candidate password against the stored Argon2 hash."""
        if not candidate or PasswordHasher is None:
            return False
        stored = self._get_password_hash()
        if not stored:
            return False
        try:
            PasswordHasher().verify(stored, candidate)
            return True
        except (VerifyMismatchError, InvalidHashError):
            return False

    # ========== Master Record System Methods ==========
    
    def _setup_master_record_cache(self) -> None:
        """Initialize cache file path when monitored_directory is known."""
        try:
            if self.monitored_directory:
                base = Path(self.monitored_directory)
                if base.exists():
                    # Store cache in the monitored directory
                    self._master_record_cache_file = base / ".samples_master_record.pkl"
                    logging.info(f"Master record cache file: {self._master_record_cache_file}")
        except Exception as e:
            logging.debug(f"Error setting up cache file: {e}")

    def _load_master_record_cache(self) -> bool:
        """Load master record from cache file. Returns True if successful."""
        try:
            if not self._master_record_cache_file or not self._master_record_cache_file.exists():
                return False
            
            with open(self._master_record_cache_file, 'rb') as f:
                cached_data = pickle.load(f)
            
            if isinstance(cached_data, dict):
                # Convert dicts back to SampleRecord objects
                with self._samples_record_lock:
                    self._samples_master_record = {
                        sid: SampleRecord(**data) if isinstance(data, dict) else data
                        for sid, data in cached_data.items()
                    }
                logging.info(f"Loaded {len(self._samples_master_record)} samples from cache")
                # Mark all as dirty to trigger UI refresh
                for record in self._samples_master_record.values():
                    record._dirty = True
                return True
        except Exception as e:
            logging.debug(f"Error loading master record cache: {e}")
        return False

    def _save_master_record_cache(self) -> None:
        """Save master record to cache file."""
        try:
            if not self._master_record_cache_file:
                return
            
            with self._samples_record_lock:
                # Convert SampleRecord objects to dicts for pickle
                cache_data = {
                    sid: asdict(record)
                    for sid, record in self._samples_master_record.items()
                }
            
            # Save to temporary file first, then rename (atomic operation)
            temp_file = self._master_record_cache_file.with_suffix('.pkl.tmp')
            with open(temp_file, 'wb') as f:
                pickle.dump(cache_data, f)
            temp_file.replace(self._master_record_cache_file)
            logging.debug(f"Saved master record cache ({len(cache_data)} samples)")
        except Exception as e:
            logging.debug(f"Error saving master record cache: {e}")

    def _ensure_job_baseline(self, sid: str) -> Dict[str, int]:
        """
        Freeze per-process baseline totals from master.csv the first time we see a sample.

        Coordinator job totals are session-only; merged lifetime totals are baseline + session.
        Baseline is captured once per process so later persisted CSV updates do not inflate it.
        """
        if sid in self._samples_job_baseline:
            return self._samples_job_baseline[sid]
        bt = bc = bf = 0
        base = self.monitored_directory
        if base:
            csv_path = Path(base) / sid / "master.csv"
            if csv_path.exists():
                try:
                    with csv_path.open("r", newline="") as fh:
                        reader = csv.DictReader(fh)
                        row = next(reader, None)
                    if row:
                        bt = int(row.get("samples_overview_total_jobs", 0) or 0)
                        bc = int(row.get("samples_overview_completed_jobs", 0) or 0)
                        bf = int(row.get("samples_overview_failed_jobs", 0) or 0)
                except Exception:
                    pass
        self._samples_job_baseline[sid] = {"total": bt, "completed": bc, "failed": bf}
        return self._samples_job_baseline[sid]

    def _merge_workflow_sample_job_counts(
        self, sid: str, workflow: Dict[str, Any]
    ) -> Dict[str, int]:
        """Merge session-only coordinator counts with frozen on-disk baseline for this process."""
        b = self._ensure_job_baseline(sid)
        wt = int(workflow.get("total_jobs", 0) or 0)
        wc = int(workflow.get("completed_jobs", 0) or 0)
        wf = int(workflow.get("failed_jobs", 0) or 0)
        wa = int(workflow.get("active_jobs", 0) or 0)
        wp = int(workflow.get("pending_jobs", 0) or 0)
        return {
            "total_jobs": b["total"] + wt,
            "completed_jobs": b["completed"] + wc,
            "failed_jobs": b["failed"] + wf,
            "active_jobs": wa,
            "pending_jobs": wp,
        }

    def _is_target_bam_finalize_redundant(self, sample_id: str) -> bool:
        """
        True when target.bam is already merged/indexed and no batch_*.bam remain.
        In that case running finalize again only queues noise jobs and confuses the samples table.
        """
        if not self.monitored_directory:
            return False
        sample_dir = Path(self.monitored_directory) / sample_id
        tb = sample_dir / "target.bam"
        bai = sample_dir / "target.bam.bai"
        if not (tb.exists() and tb.is_file() and bai.exists() and bai.is_file()):
            return False
        try:
            batch_bams = list(sample_dir.glob("batch_*.bam"))
        except Exception:
            return False
        return len(batch_bams) == 0

    def _seed_finalized_samples_from_disk(self) -> None:
        """Populate _finalized_samples from disk so restarts do not re-trigger finalization."""
        if not self.monitored_directory:
            return
        base = Path(self.monitored_directory)
        if not base.exists():
            return
        try:
            for d in base.iterdir():
                if not d.is_dir():
                    continue
                sid = d.name
                if self._is_target_bam_finalize_redundant(sid):
                    self._finalized_samples.add(sid)
        except Exception as e:
            logging.debug(f"Seed finalized samples from disk: {e}")

    def _merge_job_types_with_persisted(
        self, sid: str, workflow: Dict[str, Any]
    ) -> str:
        """Union persisted job_types (record or CSV) with workflow job types."""
        prev = ""
        try:
            with self._samples_record_lock:
                rec = self._samples_master_record.get(sid)
                if rec and rec.job_types:
                    prev = str(rec.job_types)
        except Exception:
            pass
        if not prev and self.monitored_directory:
            p = Path(self.monitored_directory) / sid / "master.csv"
            if p.exists():
                try:
                    with p.open("r", newline="") as fh:
                        reader = csv.DictReader(fh)
                        row = next(reader, None)
                    if row:
                        prev = str(row.get("samples_overview_job_types", "") or "")
                except Exception:
                    pass
        wts = workflow.get("job_types", [])
        if isinstance(wts, list):
            wjn = ",".join(sorted({str(x).strip() for x in wts if str(x).strip()}))
        else:
            wjn = str(wts or "").strip()
        parts: Set[str] = set()
        for chunk in (prev, wjn):
            for piece in (chunk or "").replace(", ", ",").split(","):
                p2 = piece.strip()
                if p2:
                    parts.add(p2)
        return ", ".join(sorted(parts))

    def _background_scan_samples(self) -> None:
        """Schedule a background scan without blocking the UI thread."""
        try:
            import asyncio
            asyncio.create_task(self._background_scan_samples_async())
        except RuntimeError:
            # Fallback if no event loop is available
            self._background_scan_samples_sync()
        except Exception as e:
            logging.debug(f"Failed to schedule background scan: {e}")

    async def _background_scan_samples_async(self) -> None:
        """Run the background scan off the UI thread and refresh UI afterwards."""
        if getattr(self, "_background_scan_in_progress", False):
            return
        self._background_scan_in_progress = True
        try:
            import asyncio
            updated_any = await asyncio.to_thread(self._background_scan_samples_sync)
            if updated_any and hasattr(self, "samples_table"):
                self._refresh_table_from_master()
        except Exception as e:
            logging.debug(f"Background scan task failed: {e}")
        finally:
            self._background_scan_in_progress = False

    def _background_scan_samples_sync(self) -> bool:
        """Scan folder and update master record (file I/O). Returns True if updated."""
        try:
            if not self.monitored_directory:
                return False
            
            base = Path(self.monitored_directory)
            if not base.exists():
                return False
            
            # Scan all sample directories
            now_ts = time.time()
            updated_any = False
            
            with self._samples_record_lock:
                # Track which samples we've seen in this scan
                seen_sample_ids = set()
                
                for sample_dir in base.iterdir():
                    if not sample_dir.is_dir():
                        continue
                    
                    sid = sample_dir.name
                    seen_sample_ids.add(sid)
                    master_csv = sample_dir / "master.csv"
                    
                    if not master_csv.exists():
                        continue
                    
                    try:
                        # Get file modification time
                        file_mtime = master_csv.stat().st_mtime
                        
                        # Get or create record
                        record = self._samples_master_record.get(sid)
                        if record is None:
                            # New sample - create record
                            record = SampleRecord(sample_id=sid)
                            record._dirty = True
                            updated_any = True
                            # Check if it's pre-existing
                            if sid in self._preexisting_sample_ids:
                                record.origin = "Pre-existing"
                            else:
                                record.origin = "Live"
                        elif file_mtime > record._file_mtime:
                            # File has changed - update record
                            record._dirty = True
                            updated_any = True
                        
                        # Update file mtime
                        record._file_mtime = file_mtime
                        
                        # Read master.csv for metadata
                        try:
                            with master_csv.open("r", newline="") as fh:
                                reader = csv.DictReader(fh)
                                first_row = next(reader, None)
                            if first_row:
                                # Update run info
                                record.run_start = self._format_timestamp_for_display(
                                    first_row.get("run_info_run_time", "")
                                )
                                record.device = first_row.get("run_info_device", "") or ""
                                record.flowcell = first_row.get("run_info_flow_cell", "") or ""
                                
                                # Update last_seen from saved value or use file mtime
                                try:
                                    saved_last = float(
                                        first_row.get("samples_overview_last_seen", 0.0) or 0.0
                                    )
                                    if saved_last > 0:
                                        record._last_seen_raw = saved_last
                                    else:
                                        record._last_seen_raw = file_mtime
                                except Exception:
                                    record._last_seen_raw = file_mtime
                                
                                # Update job counts from persisted overview
                                record.active_jobs = int(
                                    first_row.get("samples_overview_active_jobs", 0) or 0
                                )
                                record.pending_jobs = int(
                                    first_row.get("samples_overview_pending_jobs", 0) or 0
                                )
                                record.total_jobs = int(
                                    first_row.get("samples_overview_total_jobs", 0) or 0
                                )
                                record.completed_jobs = int(
                                    first_row.get("samples_overview_completed_jobs", 0) or 0
                                )
                                record.failed_jobs = int(
                                    first_row.get("samples_overview_failed_jobs", 0) or 0
                                )
                                record.job_types = str(
                                    first_row.get("samples_overview_job_types", "") or ""
                                )
                        except Exception as e:
                            logging.debug(f"Error reading master.csv for {sid}: {e}")
                        
                        # Update last_seen formatted string
                        record.last_seen = time.strftime(
                            "%Y-%m-%d %H:%M:%S", time.localtime(record._last_seen_raw)
                        )
                        # Load test_id from sample_identifier_manifest.json if present
                        record.test_id = _get_test_id_from_manifest(sample_dir)
                        
                        # Do not leave new records as default "Live" when this folder is already
                        # finished on disk — that would look like a Live→Complete transition and
                        # re-run target.bam finalization on every restart.
                        if record.origin == "Live" and (now_ts - record._last_seen_raw) >= self.completion_timeout_seconds:
                            if record.active_jobs == 0 and record.pending_jobs == 0:
                                expected = self._get_expected_completion_job_types()
                                if expected:
                                    complete_on_disk = self._expected_jobs_completed(
                                        sample_dir, expected
                                    )
                                else:
                                    complete_on_disk = self._is_target_bam_finalize_redundant(
                                        sid
                                    )
                                if not complete_on_disk:
                                    complete_on_disk = self._is_target_bam_finalize_redundant(
                                        sid
                                    )
                                if complete_on_disk:
                                    record.origin = "Complete"
                                    record._dirty = True

                        # Update origin based on inactivity timeout AND active jobs status
                        prev_origin = record.origin
                        # Only mark as Complete if timeout passed AND no active jobs
                        if record.origin == "Live" and (now_ts - record._last_seen_raw) >= self.completion_timeout_seconds:
                            if record.active_jobs == 0 and record.pending_jobs == 0:
                                record.origin = "Complete"
                                record._dirty = True
                                # Trigger finalization if transitioning from Live to Complete
                                if prev_origin == "Live" and sid not in self._finalized_samples:
                                    self._trigger_target_bam_finalization(sid)
                                    self._finalized_samples.add(sid)
                            # If there are active jobs, keep as Live even if timeout passed
                        elif record.origin == "Pre-existing" and (now_ts - record._last_seen_raw) >= self.completion_timeout_seconds:
                            # Keep as Pre-existing if it was pre-existing and still inactive
                            pass
                        elif record.origin == "Complete":
                            # Reactivate if file was modified recently OR if there are active jobs
                            if (now_ts - record._last_seen_raw) < self.completion_timeout_seconds or record.active_jobs > 0 or record.pending_jobs > 0:
                                record.origin = "Live"
                                record._dirty = True
                        
                        # Save to master record
                        self._samples_master_record[sid] = record
                        
                        # Persist overview to master.csv via MasterCSVManager
                        try:
                            manager = MasterCSVManager(str(base))
                            persist_payload = {
                                "active_jobs": int(record.active_jobs),
                                "pending_jobs": int(record.pending_jobs),
                                "total_jobs": int(record.total_jobs),
                                "completed_jobs": int(record.completed_jobs),
                                "failed_jobs": int(record.failed_jobs),
                                "job_types": record.job_types,
                                "last_seen": float(record._last_seen_raw),
                            }
                            manager.update_sample_overview(sid, persist_payload)
                        except Exception as e:
                            logging.debug(f"Error persisting overview for {sid}: {e}")
                    
                    except Exception as e:
                        logging.debug(f"Error processing sample {sid}: {e}")
                
                # Remove samples that no longer exist
                to_remove = set(self._samples_master_record.keys()) - seen_sample_ids
                if to_remove:
                    for sid in to_remove:
                        del self._samples_master_record[sid]
                    updated_any = True
                    logging.debug(f"Removed {len(to_remove)} deleted samples from master record")
            
            # Save cache if anything changed
            if updated_any:
                self._save_master_record_cache()
            
            return updated_any
        
        except Exception as e:
            logging.error(f"Error in background sample scan: {e}")
            return False

    def _refresh_table_from_master(self) -> None:
        """Fast UI refresh - reads from master record and applies filters.
        This is called from the background scanner and does NOT do file I/O."""
        try:
            if not hasattr(self, "samples_table"):
                return
            
            # Get all records from master
            with self._samples_record_lock:
                records = list(self._samples_master_record.values())
            
            # Convert to row dicts
            rows = [record.to_dict() for record in records]
            
            # Update cached rows
            self._last_samples_rows = rows
            self._last_cache_time = time.time()
            
            # Apply filters and update table
            self._apply_samples_table_filters()
            
        except Exception as e:
            logging.error(f"Error refreshing table from master: {e}")

    def _update_master_record_from_workflow(self, samples_data: List[Dict[str, Any]]) -> None:
        """Update master record from workflow polling data.
        This merges workflow stats into the master record without doing file I/O."""
        try:
            with self._samples_record_lock:
                for s in samples_data:
                    sid = s.get("sample_id", "") or "unknown"
                    if sid == "unknown":
                        continue
                    
                    last_seen = float(s.get("last_seen", time.time()))
                    
                    # Get or create record
                    record = self._samples_master_record.get(sid)
                    if record is None:
                        # Create new record (will be populated by background scan)
                        record = SampleRecord(
                            sample_id=sid,
                            _last_seen_raw=last_seen,
                            last_seen=time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(last_seen)),
                        )
                        if sid in self._preexisting_sample_ids:
                            record.origin = "Pre-existing"
                        else:
                            record.origin = "Live"
                    
                    # Merge session-only coordinator stats with on-disk baseline (see _ensure_job_baseline)
                    merged = self._merge_workflow_sample_job_counts(sid, s)
                    record.active_jobs = merged["active_jobs"]
                    record.pending_jobs = merged["pending_jobs"]
                    record.total_jobs = merged["total_jobs"]
                    record.completed_jobs = merged["completed_jobs"]
                    record.failed_jobs = merged["failed_jobs"]
                    
                    # Update file progress from job counts (same data source as other columns)
                    record.files_seen = record.total_jobs
                    record.files_processed = record.completed_jobs
                    
                    # Union with persisted job types so restarts do not drop history
                    record.job_types = self._merge_job_types_with_persisted(sid, s)
                    # Load test_id from sample_identifier_manifest.json if present
                    if self.monitored_directory:
                        sample_dir = Path(self.monitored_directory) / sid
                        record.test_id = _get_test_id_from_manifest(sample_dir)
                    
                    # Update last_seen if this is newer
                    if last_seen > record._last_seen_raw:
                        record._last_seen_raw = last_seen
                        record.last_seen = time.strftime(
                            "%Y-%m-%d %H:%M:%S", time.localtime(last_seen)
                        )
                        record._dirty = True
                    
                    # Mark as dirty to trigger UI update
                    record._dirty = True
                    self._samples_master_record[sid] = record
                
                # Persist updates to master.csv
                try:
                    base = Path(self.monitored_directory) if self.monitored_directory else None
                    if base and base.exists():
                        manager = MasterCSVManager(str(base))
                        for record in self._samples_master_record.values():
                            if record._dirty:
                                persist_payload = {
                                    "active_jobs": int(record.active_jobs),
                                    "pending_jobs": int(record.pending_jobs),
                                    "total_jobs": int(record.total_jobs),
                                    "completed_jobs": int(record.completed_jobs),
                                    "failed_jobs": int(record.failed_jobs),
                                    "job_types": record.job_types,
                                    "last_seen": float(record._last_seen_raw),
                                }
                                manager.update_sample_overview(record.sample_id, persist_payload)
                                record._dirty = False
                except Exception as e:
                    logging.debug(f"Error persisting workflow updates: {e}")
            
            # Trigger UI refresh
            if hasattr(self, "samples_table"):
                self._refresh_table_from_master()
        
        except Exception as e:
            logging.error(f"Error updating master record from workflow: {e}")

    def launch_gui(
        self,
        workflow_runner: Any = None,
        workflow_steps: list = None,
        monitored_directory: str = "",
        reload: bool = False,
        center: str = None,
    ) -> bool:
        """Launch the GUI in a completely isolated background thread.
        
        Note: When reload=True, the GUI must run in the main thread because
        signal handlers can only be set in the main thread. In this case,
        this method will block until the GUI is stopped.
        """
        if ui is None:
            logging.error("NiceGUI is not available")
            return False

        if not ensure_gui_password_set():
            logging.error("GUI password check failed. Cannot start GUI.")
            return False

        self.workflow_runner = workflow_runner
        self.workflow_steps = workflow_steps or []
        self.center = center

        # Store absolute monitored directory to avoid relative path issues
        try:
            self.monitored_directory = (
                str(Path(monitored_directory).resolve()) if monitored_directory else ""
            )
        except Exception:
            self.monitored_directory = monitored_directory
        
        # Setup master record cache when directory is known
        if self.monitored_directory:
            self._setup_master_record_cache()
            # Try to load cache immediately for fast startup
            if self._load_master_record_cache():
                logging.info("Loaded samples from cache - table will populate immediately")
            # Skip spurious target.bam finalization on restart when outputs already exist
            self._seed_finalized_samples_from_disk()

        try:
            # When reload is requested, we must run in the main thread
            # because signal handlers can only be set in the main thread
            if self.reload:
                # Check if we're in the main thread
                if threading.current_thread() is threading.main_thread():
                    # Run directly in main thread (blocking)
                    logging.info("Running GUI with reload in main thread")
                    self.is_running = True
                    try:
                        self._run_gui_worker()
                    except KeyboardInterrupt:
                        logging.info("GUI stopped by user")
                    finally:
                        self.is_running = False
                    return True
                else:
                    # Not in main thread - can't use reload
                    logging.warning(
                        "Reload requested but not in main thread. "
                        "Disabling reload functionality. "
                        "Run from main thread to enable reload."
                    )
                    self.reload = False
            
            # Start GUI in completely isolated background thread (when reload is False)
            self.gui_thread = threading.Thread(
                target=self._run_gui_worker, daemon=True, name="robin-GUI-Thread"
            )
            self.gui_thread.start()

            # NOTE: Update processing now happens inside the UI thread via ui.timer for thread-safety

            # Wait a moment for GUI to start
            time.sleep(1)

            # Check if thread is still running
            if self.gui_thread.is_alive():
                self.is_running = True
                logging.info(
                    f"GUI launched successfully on http://{self.host}:{self.port}"
                )
                return True
            else:
                logging.error("GUI thread failed to start")
                return False

        except Exception as e:
            logging.error(f"Failed to launch GUI: {e}")
            return False

    def send_update(
        self, update_type: UpdateType, data: Dict[str, Any], priority: int = 0
    ):
        """Send an update to the GUI without blocking the workflow."""
        try:
            current_time = time.time()
            queue_size = self.update_queue.qsize()
            # Adaptive threshold: higher threshold for large workloads
            # Threshold increases based on queue size to prevent dropping updates during bursts
            threshold = 500 if queue_size < 200 else 1000
            
            # Rate limiting: Skip low-priority updates if queue is getting too large
            if queue_size > threshold and priority < 5:
                # Only log when queue size changes significantly to reduce log spam
                if abs(queue_size - self._last_queue_size_logged) > 50 or (current_time - self._last_log_time) > 10.0:
                    logging.warning(
                        f"[GUI] Skipping low-priority update due to queue size: "
                        f"{queue_size} (threshold: {threshold})"
                    )
                    self._last_queue_size_logged = queue_size
                    self._last_log_time = current_time
                return
            
            # Update coalescing: Skip duplicate low-priority updates of the same type within 0.5 seconds
            # This prevents queue buildup from rapid duplicate updates
            last_time = self._last_update_times.get(update_type, 0)
            time_since_last = current_time - last_time
            
            if priority < 5 and time_since_last < 0.5:
                # Skip duplicate low-priority updates within coalescing window
                return
            
            # Create and enqueue the current update
            update = GUIUpdate(
                update_type=update_type,
                timestamp=current_time,
                data=data,
                priority=priority,
            )
            
            self._last_update_times[update_type] = current_time

            # Use negative priority so higher priority updates come first.
            # Include a monotonically increasing sequence as a tiebreaker so
            # heap comparisons never fallback to comparing GUIUpdate objects.
            self._update_seq += 1
            self.update_queue.put((-priority, self._update_seq, update))
            self.total_updates_enqueued += 1
            
            # Reduced logging: only log every 100th update or important updates
            if self.total_updates_enqueued % 100 == 0 or priority >= 8:
                logging.debug(
                    f"[GUI] Enqueued update #{self.total_updates_enqueued}: {update_type.value} (queue: {queue_size})"
                )

        except Exception as e:
            # Don't let update failures affect the workflow
            logging.debug(f"Failed to send GUI update: {e}")

    def _process_progress_queue(self):
        """Process queued progress updates for report generation."""
        try:
            from robin.gui.report_progress import progress_manager
            
            # Process updates in the UI context
            while not progress_manager.progress_queue.empty():
                update = progress_manager.progress_queue.get_nowait()
                self._handle_progress_update(update)
                
        except queue.Empty:
            pass
        except Exception as e:
            logging.debug(f"Error processing progress queue: {e}")
    
    def _setup_notification_system(self, container):
        """Set up the notification system with the provided container."""
        try:
            # Set up event subscription in the proper UI context
            self.progress_notification_event.subscribe(
                lambda data: self._show_notification_in_container(container, data)
            )
            logging.info("Notification system set up successfully")
        except Exception as e:
            logging.error(f"Error setting up notification system: {e}")

    def _show_notification_in_container(self, container, data: Dict[str, Any]):
        """Show notification in the dedicated container."""
        try:
            message = data.get('message', '')
            notification_type = data.get('type', 'info')
            timeout = data.get('timeout', 5000)
            
            # Create notification in the container
            with container:
                ui.notify(
                    message,
                    type=notification_type,
                    timeout=timeout,
                    position="top-right"
                )
                
        except Exception as e:
            logging.error(f"Error showing notification in container: {e}")

    def _handle_progress_notification_event(self, event_data: Dict[str, Any]):
        """Handle progress notification events."""
        try:
            from nicegui import ui

            message = event_data.get('message', '')
            notification_type = event_data.get('type', 'info')
            timeout = event_data.get('timeout', 5000)

            # Show notification in the proper UI context
            ui.notify(
                message,
                type=notification_type,
                timeout=timeout,
                position="top-right"
            )

        except Exception as e:
            logging.error(f"Error handling progress notification event: {e}")

    def _handle_progress_update(self, update: Dict[str, Any]):
        """Handle a single progress update in the UI context."""
        try:
            from nicegui import ui
            
            sample_id = update['sample_id']
            update_type = update['type']
            
            if update_type == 'start':
                # Emit event for initial notification
                self.progress_notification_event.emit({
                    'sample_id': sample_id,
                    'message': f"[{sample_id}] Starting report generation...",
                    'type': 'info',
                    'timeout': 3000  # 3 seconds
                })
                
                logging.debug(f"Started report generation for {sample_id}")
                
            elif update_type == 'update':
                stage = update['stage']
                message = update['message']
                progress = update.get('progress')
                
                # Calculate progress percentage
                progress_percent = int((progress or 0.0) * 100) if progress is not None else ""
                progress_text = f" ({progress_percent}%)" if progress_percent else ""
                
                # Emit event for progress notification
                self.progress_notification_event.emit({
                    'sample_id': sample_id,
                    'message': f"[{sample_id}] {message}{progress_text}",
                    'type': 'info',
                    'timeout': 2000  # 2 seconds for progress updates
                })
                
                logging.debug(f"Updated progress for {sample_id}: {stage} - {message}")
                
            elif update_type == 'complete':
                filename = update.get('filename')
                
                # Show completion notification
                completion_message = f"[{sample_id}] Report generation completed"
                if filename:
                    completion_message += f": {filename}"
                
                # Emit event for completion notification
                self.progress_notification_event.emit({
                    'sample_id': sample_id,
                    'message': completion_message,
                    'type': 'positive',
                    'timeout': 5000
                })
                
                logging.info(f"Completed report generation for {sample_id}")
                
            elif update_type == 'error':
                error_message = update['error_message']
                
                # Emit event for error notification
                self.progress_notification_event.emit({
                    'sample_id': sample_id,
                    'message': f"[{sample_id}] Report generation failed: {error_message}",
                    'type': 'negative',
                    'timeout': 10000
                })
                
                logging.error(f"Report generation failed for {sample_id}: {error_message}")
                
        except Exception as e:
            logging.error(f"Error handling progress update: {e}")

    def _drain_updates_on_ui(self):
        """Drain queued updates and apply them on the UI thread (called by ui.timer).
        
        Uses adaptive throttling: processes more updates when queue is large,
        fewer when queue is small. Batches UI updates for better performance.
        """
        if not self.gui_ready.is_set():
            return
        
        processed_count = 0
        max_updates_per_cycle = 50  # Limit updates per cycle to prevent UI blocking
        
        try:
            queue_size = self.update_queue.qsize()
            
            # Adaptive processing: process more updates when queue is large
            if queue_size > 200:
                max_updates_per_cycle = 100  # Process more aggressively when backlogged
            elif queue_size > 100:
                max_updates_per_cycle = 75
            elif queue_size < 20:
                max_updates_per_cycle = 25  # Process fewer when queue is small
            
            # Process updates in batch
            while processed_count < max_updates_per_cycle:
                try:
                    # Entries are (-priority, seq, GUIUpdate)
                    _, _, update = self.update_queue.get_nowait()
                except queue.Empty:
                    break
                
                self._handle_update(update)
                self.total_updates_processed += 1
                processed_count += 1
                
                # Reduced logging: only log every 50th update or when queue size changes significantly
                if self.total_updates_processed % 50 == 0:
                    logging.debug(
                        f"[GUI] Processed update #{self.total_updates_processed}: {update.update_type.value} (queue: {queue_size})"
                    )
            
            # Note: Timer interval is fixed at 0.5s for consistent performance
            # Adaptive processing (max_updates_per_cycle) handles queue size variations
                
        except Exception as e:
            logging.debug(f"[GUI] Error draining updates on UI: {e}")
        finally:
            # Batch UI update: only call ui.update() once per drain cycle
            if processed_count > 0:
                try:
                    ui.update()
                except Exception:
                    pass

    def _handle_update(self, update: GUIUpdate):
        """Handle a single update message."""
        try:
            if update.update_type == UpdateType.WORKFLOW_STATUS:
                self._update_workflow_status(update.data)
            elif update.update_type == UpdateType.JOB_UPDATE:
                self._update_job_status(update.data)
            elif update.update_type == UpdateType.QUEUE_UPDATE:
                self._update_queue_status(update.data)
            elif update.update_type == UpdateType.LOG_MESSAGE:
                self._update_logs(update.data)
            elif update.update_type == UpdateType.PROGRESS_UPDATE:
                self._update_progress(update.data)
            elif update.update_type == UpdateType.ERROR_UPDATE:
                self._update_errors(update.data)
            elif update.update_type == UpdateType.SAMPLES_UPDATE:
                # Update master record from workflow polling data
                # This merges workflow stats into master record without file I/O
                self._pending_samples_data = update.data
                try:
                    samples_list = update.data.get("samples", [])
                    if samples_list:
                        self._update_master_record_from_workflow(samples_list)
                except Exception as e:
                    logging.error(f"[GUI] Samples master record update failed: {e}")
            elif update.update_type == UpdateType.WARNING_NOTIFICATION:
                self._show_warning_notification(update.data)

        except Exception as e:
            logging.debug(f"Error handling GUI update: {e}")

    def _update_workflow_status(self, data: Dict[str, Any]):
        """Update workflow status in the GUI."""
        try:
            if hasattr(self, "status_indicator") and hasattr(self, "status_label"):
                if data.get("is_running", False):
                    self.status_indicator.classes(
                        replace=(
                            "text-2xl workflow-monitor-status-icon "
                            "workflow-monitor-status-icon--running"
                        )
                    )
                    self.status_label.set_text("Workflow status: Running")
                    self.status_label.classes(
                        replace="text-sm font-medium workflow-monitor-status-text--running"
                    )
                    self._is_running = True
                else:
                    self.status_indicator.classes(
                        replace=(
                            "text-2xl workflow-monitor-status-icon "
                            "workflow-monitor-status-icon--stopped"
                        )
                    )
                    self.status_label.set_text("Workflow status: Stopped")
                    self.status_label.classes(
                        replace="text-sm font-medium workflow-monitor-status-text--stopped"
                    )
                    self._is_running = False

                # Update timing
                if data.get("start_time"):
                    self._start_time = float(data["start_time"])
                    if hasattr(self, "workflow_start_time"):
                        start_str = time.strftime(
                            "%Y-%m-%d %H:%M:%S", time.localtime(self._start_time)
                        )
                        self.workflow_start_time.set_text(f"Started: {start_str}")

                if hasattr(self, "workflow_duration") and self._start_time:
                    elapsed_seconds = int(time.time() - self._start_time)
                    self.workflow_duration.set_text(
                        f"Duration: {self._format_duration(elapsed_seconds)}"
                    )

        except Exception as e:
            logging.debug(f"Error updating workflow status: {e}")

    def _update_job_status(self, data: Dict[str, Any]):
        """Update job status in the GUI with search/sort/filter support."""
        try:
            if hasattr(self, "active_jobs_table"):
                # Build all rows
                job_rows: List[Dict[str, Any]] = []
                for job in data.get("active_jobs", []):
                    job_rows.append(
                        {
                            "job_id": str(job.get("job_id", ""))[:8],
                            "job_type": job.get("job_type", ""),
                            "filepath": Path(job.get("filepath", "")).name,
                            "worker": job.get("worker_name", ""),
                            "duration": self._format_duration(job.get("duration", 0)),
                            "progress": f"{int(job.get('progress', 0) * 100)}%",
                        }
                    )

                # Cache for filtering
                self._active_jobs_all_rows = job_rows

                # Update filter options dynamically
                try:
                    type_options = ["All"] + sorted(
                        {r.get("job_type", "") for r in job_rows if r.get("job_type")}
                    )
                    worker_options = ["All"] + sorted(
                        {r.get("worker", "") for r in job_rows if r.get("worker")}
                    )
                    if hasattr(self, "active_jobs_type_filter"):
                        self.active_jobs_type_filter.set_options(type_options)
                    if hasattr(self, "active_jobs_worker_filter"):
                        self.active_jobs_worker_filter.set_options(worker_options)
                except Exception:
                    pass

                # Apply filters and refresh table
                self._apply_active_jobs_filters_and_update()

                # Show/hide the "no active jobs" placeholder
                if hasattr(self, "no_active_jobs_label"):
                    if job_rows:
                        self.no_active_jobs_label.set_visibility(False)
                    else:
                        self.no_active_jobs_label.set_visibility(True)

        except Exception as e:
            logging.debug(f"Error updating job status: {e}")

    def _update_queue_status(self, data: Dict[str, Any]):
        """Update queue status in the GUI."""
        try:
            # Cache the latest payload so newly created pages can show current values
            self._last_queue_status = data or {}
            # Emit a concise log line to the Live Logs pane for visibility
            try:
                pr = data.get("preprocessing", {})
                ar = {
                    "mgmt": data.get("mgmt", {}),
                    "cnv": data.get("cnv", {}),
                    "target": data.get("target", {}),
                    "fusion": data.get("fusion", {}),
                }
                an_running = sum(int(v.get("running", 0) or 0) for v in ar.values())
                an_total = sum(int(v.get("total", 0) or 0) for v in ar.values())
                cl = data.get("classification", {})
                ot = data.get("other", {})
                msg = (
                    f"Queues | Pre:{pr.get('running',0)}/{pr.get('total',0)} "
                    f"| An:{an_running}/{an_total} "
                    f"| Cl:{cl.get('running',0)}/{cl.get('total',0)} "
                    f"| Ot:{ot.get('running',0)}/{ot.get('total',0)}"
                )
                self._log_buffer.append(f"[queue] {msg}\n")
                if hasattr(self, "log_area"):
                    self.log_area.set_value("".join(self._log_buffer))
            except Exception:
                pass
            # Update queue status displays
            if hasattr(self, "preprocessing_status"):
                queue_data = data.get("preprocessing", {})
                self.preprocessing_status.set_text(
                    f"{queue_data.get('running', 0)}/{queue_data.get('total', 0)}"
                )

            if hasattr(self, "analysis_status"):
                # Combine analysis queues
                analysis_running = sum(
                    data.get(q, {}).get("running", 0)
                    for q in ["mgmt", "cnv", "target", "fusion"]
                )
                analysis_total = sum(
                    data.get(q, {}).get("total", 0)
                    for q in ["mgmt", "cnv", "target", "fusion"]
                )
                self.analysis_status.set_text(f"{analysis_running}/{analysis_total}")

            if hasattr(self, "classification_status"):
                classification_data = data.get("classification", {})
                self.classification_status.set_text(
                    f"{classification_data.get('running', 0)}/{classification_data.get('total', 0)}"
                )

            if hasattr(self, "other_status"):
                other_data = data.get("other", {})
                self.other_status.set_text(
                    f"{other_data.get('running', 0)}/{other_data.get('total', 0)}"
                )

        except Exception as e:
            logging.debug(f"Error updating queue status: {e}")

    def _apply_active_jobs_filters_and_update(self) -> None:
        """Apply dropdown filters and current search to the Active Jobs table."""
        try:
            all_rows = getattr(self, "_active_jobs_all_rows", [])
            if not hasattr(self, "active_jobs_table"):
                return
            selected_type = "All"
            selected_worker = "All"
            try:
                if hasattr(self, "active_jobs_type_filter"):
                    selected_type = self.active_jobs_type_filter.value or "All"
                if hasattr(self, "active_jobs_worker_filter"):
                    selected_worker = self.active_jobs_worker_filter.value or "All"
            except Exception:
                pass

            def _keep(row: Dict[str, Any]) -> bool:
                if selected_type != "All" and row.get("job_type") != selected_type:
                    return False
                if selected_worker != "All" and row.get("worker") != selected_worker:
                    return False
                return True

            filtered = [r for r in all_rows if _keep(r)]
            self.active_jobs_table.rows = filtered
            self.active_jobs_table.update()
        except Exception:
            pass

    def _update_logs(self, data: Dict[str, Any]):
        """Update logs in the GUI."""
        try:
            if hasattr(self, "log_area"):
                log_message = data.get("message", "")
                log_level = data.get("level", "INFO")
                timestamp = time.strftime("%H:%M:%S")

                new_line = f"[{timestamp}] {log_level}: {log_message}\n"
                self._log_buffer.append(new_line)
                self.log_area.set_value("".join(self._log_buffer))

        except Exception as e:
            logging.debug(f"Error updating logs: {e}")

    def _update_progress(self, data: Dict[str, Any]):
        """Update progress in the GUI."""
        try:
            if hasattr(self, "progress_bar") and hasattr(self, "progress_label"):
                progress = data.get("progress", 0.0)

                pct = max(0.0, min(100.0, round(progress * 100.0, 1)))

                self.progress_bar.set_value(round(progress, 2))
                # Drop .0 for integers like 81.0 -> 81
                pct_str = f"{pct:.1f}" if pct % 1 else f"{int(pct)}"
                self.progress_label.set_text(f"{pct_str}% Complete")
                # Optional counts
                if hasattr(self, "completed_count") and "completed" in data:
                    self.completed_count.set_text(str(data["completed"]))
                if hasattr(self, "failed_count") and "failed" in data:
                    self.failed_count.set_text(str(data["failed"]))
                if hasattr(self, "total_count") and "total" in data:
                    self.total_count.set_text(str(data["total"]))

        except Exception as e:
            logging.debug(f"Error updating progress: {e}")

    def _update_file_progress(self):
        """Update file progress display in workflow monitor from current samples data."""
        try:
            if not hasattr(self, "sample_files_progress_container") or not hasattr(self, "_last_samples_rows"):
                return
            
            rows = self._last_samples_rows or []
            total_files_seen = 0
            total_files_processed = 0
            sample_progress_data = []
            
            # Calculate totals and collect per-sample data
            for row in rows:
                files_seen = row.get("files_seen", 0) or 0
                files_processed = row.get("files_processed", 0) or 0
                sample_id = row.get("sample_id", "")
                
                if files_seen > 0:
                    total_files_seen += files_seen
                    total_files_processed += files_processed
                    sample_progress_data.append({
                        "sample_id": sample_id,
                        "files_seen": files_seen,
                        "files_processed": files_processed,
                        "progress": files_processed / files_seen if files_seen > 0 else 0.0
                    })
            
            # Update overall progress
            if hasattr(self, "overall_files_progress") and hasattr(self, "overall_files_label"):
                overall_progress = total_files_processed / total_files_seen if total_files_seen > 0 else 0.0
                self.overall_files_progress.set_value(round(overall_progress, 2))
                self.overall_files_label.set_text(f"{total_files_processed}/{total_files_seen} files processed")
            
            # Update per-sample progress bars (limit to first 20 to avoid UI overload)
            if hasattr(self, "sample_files_progress_container"):
                try:
                    # Clear existing content
                    self.sample_files_progress_container.clear()
                    
                    if sample_progress_data:
                        # Sort by progress (lowest first) to show samples needing attention
                        sample_progress_data.sort(key=lambda x: x["progress"])
                        
                        # Show up to 20 samples with most activity
                        for sample_info in sample_progress_data[:20]:
                            sample_id = sample_info["sample_id"]
                            progress = sample_info["progress"]
                            files_seen = sample_info["files_seen"]
                            files_processed = sample_info["files_processed"]
                            
                            with self.sample_files_progress_container:
                                with ui.row().classes(
                                    "w-full items-center gap-2 mb-1 min-w-0"
                                ):
                                    ui.label(sample_id).classes(
                                        "text-xs font-mono min-w-[120px] "
                                        "workflow-monitor-sample-id"
                                    )
                                    progress_bar = ui.linear_progress(progress).classes(
                                        "flex-1 min-w-0"
                                    )
                                    if progress >= 1.0:
                                        progress_bar.props("color=positive")
                                    ui.label(f"{files_processed}/{files_seen}").classes(
                                        "text-xs min-w-[60px] workflow-monitor-file-count"
                                    )
                    else:
                        # Show placeholder when no data
                        with self.sample_files_progress_container:
                            ui.label(
                                "No file progress data available yet."
                            ).classes("classification-insight-foot italic")
                except Exception as e:
                    logging.debug(f"Error updating per-sample file progress: {e}")
                            
        except Exception as e:
            logging.debug(f"Error updating file progress: {e}")

    def _update_errors(self, data: Dict[str, Any]):
        """Update error information in the GUI."""
        try:
            # Update error counts if available
            if hasattr(self, "preprocessing_errors"):
                self.preprocessing_errors.set_text(
                    str(data.get("preprocessing_errors", 0))
                )

            if hasattr(self, "analysis_errors"):
                self.analysis_errors.set_text(str(data.get("analysis_errors", 0)))

            if hasattr(self, "classification_errors"):
                self.classification_errors.set_text(
                    str(data.get("classification_errors", 0))
                )

        except Exception as e:
            logging.debug(f"Error updating error information: {e}")

    def _show_warning_notification(self, data: Dict[str, Any]):
        """Show a warning notification in the GUI with dismiss button."""
        try:
            if ui is None:
                return
            
            message = data.get("message", "")
            title = data.get("title", "Warning")
            sample_id = data.get("sample_id", "")
            filename = data.get("filename", "")
            level = data.get("level", "warning")
            
            # Build notification message
            if sample_id:
                notification_msg = f"[{sample_id}] {message}"
            elif filename:
                notification_msg = f"[{filename}] {message}"
            else:
                notification_msg = message
            
            # Show dismissible notification
            # NiceGUI notifications are dismissible by default with a close button
            # timeout=0 makes it persistent until manually dismissed
            ui.notify(
                notification_msg,
                type=level,
                timeout=0,  # Persistent until manually dismissed (close button available)
                position="top-right",
            )
            
        except Exception as e:
            logging.debug(f"Error showing warning notification: {e}")

    def _format_duration(self, seconds):
        """Format duration in seconds to human readable format."""
        if seconds < 60:
            return f"{seconds}s"
        elif seconds < 3600:
            minutes = seconds // 60
            return f"{minutes}m {seconds % 60}s"
        else:
            hours = seconds // 3600
            minutes = (seconds % 3600) // 60
            return f"{hours}h {minutes}m"

    def _run_gui_worker(self):
        """Run the GUI in a completely isolated thread."""
        try:
            # Set thread name for identification
            threading.current_thread().name = "robin-GUI-Thread"
            self._register_auth_middleware()

            # Session cookie dark_mode (nicegui.io pattern): available on first HTML response so theme matches paint.
            if app is not None:

                @app.post("/robin_dark_mode")
                async def _robin_dark_mode(request: Request):
                    data = await request.json()
                    v = data.get("value")
                    app.storage.browser["dark_mode"] = bool(v)
                    try:
                        app.storage.user["dark_mode"] = bool(v)
                    except Exception:
                        pass

            @ui.page("/login")
            def login_page(redirect_to: str = "/"):
                """Authenticate user session before allowing access."""
                if self._is_authenticated():
                    safe_target = redirect_to if redirect_to and redirect_to != "/login" else "/"
                    return RedirectResponse(safe_target)

                def try_login() -> None:
                    if self._verify_password(password.value):
                        app.storage.user.update({
                            "authenticated": True,
                            "_auth_generation": app.storage.general.get("_auth_generation"),
                        })
                        safe_target = (
                            redirect_to if redirect_to and redirect_to != "/login" else "/"
                        )
                        ui.navigate.to(safe_target)
                    else:
                        ui.notify("Incorrect password", type="negative")

                _setup_global_resources()
                with theme.frame(
                    "<strong>R</strong>apid nanop<strong>O</strong>re <strong>B</strong>rain intraoperat<strong>I</strong>ve classificatio<strong>N</strong>",
                    smalltitle="<strong>R.O.B.I.N</strong>",
                    batphone=False,
                    center=self.center,
                ):
                    with ui.element("div").classes("w-full min-w-0").props(
                        "id=login-page"
                    ):
                        with ui.column().classes(
                            "w-full min-h-[70vh] items-center justify-center "
                            "p-4 md:p-6"
                        ):
                            with ui.element("div").classes(
                                "w-full max-w-md classification-insight-shell "
                                "min-w-0"
                            ):
                                ui.label("Sign in").classes(
                                    "classification-insight-heading "
                                    "text-headline-small text-center w-full"
                                )
                                with ui.element("div").classes(
                                    "classification-insight-card w-full min-w-0"
                                ):
                                    with ui.column().classes(
                                        "w-full min-w-0 gap-3 p-2 md:p-3"
                                    ):
                                        with ui.row().classes(
                                            "items-center gap-2 justify-center "
                                            "flex-nowrap min-w-0"
                                        ):
                                            ui.icon("lock").classes(
                                                "classification-insight-icon shrink-0"
                                            )
                                            ui.label("R.O.B.I.N").classes(
                                                "classification-insight-model "
                                                "shrink-0 text-center"
                                            )
                                        ui.label(
                                            "Enter password to continue."
                                        ).classes(
                                            "classification-insight-foot text-center"
                                        )
                                        password = (
                                            ui.input(
                                                "Password",
                                            )
                                            .props(
                                                "autocomplete=off outlined dense"
                                            )
                                            .classes("w-full")
                                            .style("-webkit-text-security: disc;")
                                            .on("keydown.enter", try_login)
                                        )
                                        ui.button(
                                            "Log in",
                                            on_click=try_login,
                                            icon="login",
                                        ).props("color=primary no-caps").classes(
                                            "w-full"
                                        )

            # Create the main workflow monitor page
            @ui.page("/")
            def welcome_page():
                """Welcome page at root route."""
                _setup_global_resources()
                self._create_welcome_page()

            # Create the workflow monitoring page
            @ui.page("/robin")
            def workflow_monitor():
                """Workflow monitoring page under /robin route."""
                _setup_global_resources()
                self._create_workflow_monitor()

            # Create the samples overview page
            @ui.page("/live_data")
            def samples_overview():
                """Samples overview page showing all tracked samples."""
                _setup_global_resources()
                logging.info("[samples_overview] building page /live_data")
                self._create_samples_overview()

            # Create individual sample detail pages
            @ui.page("/live_data/{sample_id}")
            def sample_detail(sample_id: str):
                """Individual sample detail page."""
                _setup_global_resources()

                # Clear cached state BEFORE creating the page so plots refresh on load
                # This ensures all graphs load properly on each page visit
                self._refresh_sample_plots(sample_id)
                
                # Add a page visit handler to refresh plots on reconnect
                def on_page_visit():
                    """Handle page visit - refresh plots if needed."""
                    try:
                        logging.info(
                            f"Page visit detected for sample {sample_id} - triggering plot refresh"
                        )
                        # Small delay to ensure components are loaded
                        ui.timer(
                            1.0,
                            lambda: self._refresh_sample_plots(sample_id),
                            once=True,
                        )
                    except Exception as e:
                        logging.debug(f"Page visit handler failed for {sample_id}: {e}")

                # Register the page visit handler
                ui.context.client.on_connect(on_page_visit)

                self._create_sample_detail_page(sample_id)

            # Create sample details page
            @ui.page("/live_data/{sample_id}/details")
            def sample_details(sample_id: str):
                """Sample details page with comprehensive information."""
                _setup_global_resources()
                self._create_sample_details_page(sample_id)

            # Watched folders management page
            @ui.page("/watched_folders")
            def watched_folders_page():
                """Page for managing watched folders."""
                _setup_global_resources()
                self._create_watched_folders_page()

            # Sample identifier generator page
            @ui.page("/sample_id_generator")
            def sample_id_generator_page():
                """Page for generating sample identifiers from Test ID, name, and D.O.B."""
                _setup_global_resources()
                self._create_sample_id_generator_page()

            # Download API endpoint
            @ui.page("/api/download/{sample_id}/{filename}")
            def download_file(sample_id: str, filename: str):
                """Download a file from a sample directory."""
                try:
                    # Security: Only allow alphanumeric characters and common file extensions
                    import re
                    if not re.match(r'^[a-zA-Z0-9._-]+$', filename):
                        ui.notify("Invalid filename", type="error")
                        return
                    
                    # Find the sample directory
                    base_dir = Path(self.monitored_directory) if self.monitored_directory else None
                    if not base_dir or not base_dir.exists():
                        ui.notify("Sample directory not found", type="error")
                        return
                    
                    sample_dir = base_dir / sample_id
                    if not sample_dir.exists():
                        ui.notify(f"Sample {sample_id} not found", type="error")
                        return
                    
                    file_path = sample_dir / filename
                    if not file_path.exists() or not file_path.is_file():
                        ui.notify(f"File {filename} not found", type="error")
                        return
                    
                    # Read file content
                    with open(file_path, 'rb') as f:
                        content = f.read()
                    
                    # Use NiceGUI's download functionality
                    ui.download(
                        content,
                        filename=filename,
                        media_type='application/octet-stream'
                    )
                    
                except Exception as e:
                    ui.notify(f"Download failed: {e}", type="error")

            # Setup global CSS and static files - moved to a helper function
            def _setup_global_resources():
                """Setup global CSS and static file resources."""
                ui.add_css(
                    """
                    .shadows-into light-regular {
                        font-family: "Shadows Into Light", cursive;
                        font-weight: 800;
                        font-style: normal;
                    }
                """
                )
                # Register fonts from the GUI package if available
                try:
                    fonts_dir = Path(__file__).parent / "gui" / "fonts"
                    if fonts_dir.exists():
                        app.add_static_files("/fonts", str(fonts_dir))
                    else:
                        logging.debug(f"Fonts directory not found: {fonts_dir}")
                except Exception as e:
                    logging.debug(f"Could not register fonts static dir: {e}")

            # Setup global timers and processing - moved to a helper function
            def _setup_global_timers():
                """Setup global timers and update processing."""
                try:
                    self.gui_ready.set()
                    # Faster drain rate: 0.5s to keep up with high-volume updates
                    # The _drain_updates_on_ui method adaptively processes more/fewer updates per cycle
                    app.timer(0.5, self._drain_updates_on_ui, active=True)
                    
                    # Background master record scanner - scans folder periodically
                    # Runs every 10 seconds to keep master record up to date
                    # First scan happens immediately if cache was loaded, otherwise after 1 second
                    initial_delay = 1.0 if not self._samples_master_record else 0.1
                    app.timer(
                        initial_delay,
                        lambda: self._background_scan_samples(),
                        once=True,
                    )
                    # Subsequent scans every 10 seconds
                    app.timer(self._background_scan_interval, self._background_scan_samples, active=True)
                    
                    # If cache was not loaded, do initial preexisting scan after GUI is ready
                    if not self._samples_master_record:
                        app.timer(
                            2.0,
                            lambda: self._background_scan_samples(),  # Initial scan
                            once=True,
                        )
                    
                    # Process progress queue for report generation
                    app.timer(0.1, self._process_progress_queue, active=True)
                    
                except Exception:
                    pass

            try:
                iconfile = os.path.join(
                    os.path.dirname(os.path.abspath(images.__file__)), "favicon.ico"
                )
                if not os.path.exists(iconfile):
                    logging.warning(f"Favicon file not found: {iconfile}")
            except Exception as e:
                logging.error(f"Error locating favicon: {str(e)}")
                iconfile = None
            
            # Setup global timers once when the app starts
            _setup_global_timers()
            
            try:
                from robin_native_api import attach_native_api
                attach_native_api(app, self)
            except Exception as e:
                logging.warning("Native API not mounted: %s", e)
            
            
            
            # Start the GUI
            ui.run(
                host=self.host,
                port=self.port,
                show=False,
                reload=self.reload,
                title="ROBIN",
                storage_secret="robin",
                favicon=iconfile,
                reconnect_timeout=60,
            )
        except Exception as e:
            print(f"GUI worker error: {e}")
            import traceback

            traceback.print_exc()

    def _create_welcome_page(self):
        """Create the welcome page (Editorial Bioinformatics, design.md)."""
        # Surfaces follow global body gradient + .q-card (slate border, light elevation).
        with theme.frame(
            "<strong>R</strong>apid nanop<strong>O</strong>re <strong>B</strong>rain intraoperat<strong>I</strong>ve classificatio<strong>N</strong>",
            smalltitle="<strong>R.O.B.I.N</strong>",
            batphone=False,
            center=self.center,
            setup_notifications=self._setup_notification_system,
        ):
            with ui.column().classes("w-full min-h-screen p-3 md:p-6"):
                # Anchor: Manrope display, slate hierarchy (design.md typography + light surfaces)
                with ui.card().classes("w-full max-w-6xl mx-auto mb-8"):
                    with ui.column().classes(
                        "p-6 gap-3 w-full bg-gradient-to-b from-slate-50 to-white "
                        "dark:from-slate-900/80 dark:to-[var(--md-surface)] rounded-[inherit]"
                    ):
                        with ui.row().classes("w-full justify-end"):
                            ui.image(theme.get_imagefile()).classes(
                                "flex-shrink-0 object-contain w-28 sm:w-32 md:w-36 "
                                "max-h-36"
                            )
                        ui.label("Welcome to R.O.B.I.N").classes(
                            "text-display-small text-slate-900 dark:text-slate-50 "
                            "text-left w-full"
                        )
                        ui.label(
                            "This tool enables real time analysis of data from Oxford Nanopore Technologies sequencers."
                        ).classes(
                            "text-body-large text-slate-600 dark:text-slate-400 "
                            "text-left w-full"
                        )

                # Modular content card: Inter body, emerald primary actions
                with ui.card().classes("w-full max-w-6xl mx-auto mb-8"):
                    with ui.column().classes("p-6 gap-4"):
                        ui.label("What is R.O.B.I.N?").classes(
                            "text-headline-medium text-center text-slate-900 dark:text-slate-50"
                        )
                        ui.label(
                            "ROBIN (Rapid nanopOre Brain intraoperatIve classificatioN) is a comprehensive bioinformatics workflow system designed for processing and analyzing BAM files. It provides automated preprocessing, multiple analysis pipelines, and real-time monitoring capabilities. It now incorporates LITTLE JOHN (Lightweight Infrastructure for Task Tracking and Logging with Extensible Job Orchestration for High-throughput aNalysis), which handles the heavy lifting behind the scenes."
                        ).classes(
                            "text-body-large text-slate-600 dark:text-slate-400 text-center"
                        )

                        _cta_primary = (
                            "bg-primary text-white rounded-lg px-6 py-4 text-title-medium "
                            "text-center transition-opacity hover:opacity-90 no-underline"
                        )
                        _cta_secondary = (
                            "border border-primary text-primary rounded-lg px-6 py-4 "
                            "text-title-medium text-center bg-transparent "
                            "hover:bg-[var(--md-primary-container)] "
                            "hover:text-[var(--md-on-primary-container)] "
                            "transition-colors no-underline"
                        )

                        with ui.grid(columns=2).classes("w-full gap-3"):
                            ui.link("View All Samples", "/live_data").classes(
                                _cta_primary
                            )
                            ui.link("Open Workflow Monitor", "/robin").classes(
                                _cta_primary
                            )
                            ui.link(
                                "Manage watched folders",
                                "/watched_folders",
                            ).classes(_cta_secondary)
                            ui.link(
                                "Generate Sample ID",
                                "/sample_id_generator",
                            ).classes(_cta_secondary)
                            ui.link(
                                "View Documentation",
                                "https://looselab.github.io/ROBIN/",
                            ).classes(_cta_primary)

                # News — card uses same surface treatment as global theme
                with ui.card().classes("w-full max-w-6xl mx-auto mb-8"):
                    with ui.column().classes("w-full"):
                        # Initialize news feed only if it hasn't been initialized yet
                        if self.news_feed is None:
                            self.news_feed = NewsFeed()
                            self.news_feed.start_update_timer()
                        # Create the news element
                        self.news_feed.create_news_element()

    def _create_samples_overview(self):
        """Create the samples overview page showing all tracked samples (design.md Editorial Bioinformatics)."""
        logging.info("[samples_overview] _create_samples_overview() started")
        with theme.frame(
            "R.O.B.I.N - Sample Tracking Overview",
            smalltitle="Samples",
            batphone=False,
            center=self.center,
            setup_notifications=self._setup_notification_system,
        ):
            with ui.column().classes("w-full min-h-[70vh] p-3 md:p-6 gap-4"):
                # Anchor: bold headline + secondary copy (Digital Curator hierarchy)
                # w-full on inner column so the gradient spans the full card (not content width)
                with ui.card().classes("w-full max-w-7xl mx-auto overflow-hidden"):
                    with ui.column().classes(
                        "w-full min-w-0 p-4 md:p-6 gap-3 bg-gradient-to-b from-slate-50 to-white "
                        "dark:from-slate-900/80 dark:to-[var(--md-surface)]"
                    ):
                        ui.label("Sample tracking").classes(
                            "text-display-small text-slate-900 dark:text-slate-50"
                        )
                        ui.label(
                            "Tracked nanopore runs and workflow state in one place. "
                            "Use search and origin filters to find a library, then open a sample for analysis."
                        ).classes(
                            "text-body-large text-slate-600 dark:text-slate-400 max-w-3xl"
                        )

                # Samples table section — outer column keeps mobile scroll behavior
                with ui.column().classes("w-full max-w-7xl mx-auto gap-3"):
                    # Title + export/SNP actions on one row (must stay visible; a separate row below filters was easy to miss)
                    with ui.row().classes(
                        "w-full items-center justify-between flex-wrap gap-3 mb-2"
                    ):
                        ui.label("All tracked samples").classes(
                            "text-headline-medium text-slate-900 dark:text-slate-50 shrink-0"
                        )
                        ui.label("A=Active  P=Pending  T=Total  C=Completed  F=Failed").classes(
                            "text-[11px] text-slate-600 dark:text-slate-400 shrink-0"
                        )
                        with ui.row().classes(
                            "gap-2 flex-wrap items-center justify-end min-w-0 flex-1"
                        ):
                            self.samples_loading_indicator = ui.spinner(
                                size="sm", color="primary"
                            )
                            self.samples_loading_indicator.set_visibility(False)
                            ui.button(
                                "Select all",
                                on_click=lambda: self._samples_select_all_visible_for_export(),
                            ).props("flat dense no-caps outline")
                            ui.button(
                                "Clear selection",
                                on_click=lambda: self._samples_clear_export_selection(),
                            ).props("flat dense no-caps outline")
                            self.bulk_snp_button = ui.button(
                                "SNP: all missing",
                                icon="biotech",
                                on_click=lambda: None,
                            ).props("color=secondary dense no-caps").classes(
                                "border border-slate-300 dark:border-slate-600"
                            )
                            # Optional bulk MNP-Flex action (only if credentials are available)
                            _mnpflex_username = (
                                os.getenv("MNPFLEX_USERNAME") or os.getenv("EPIGNOSTIX_USERNAME")
                            )
                            _mnpflex_password = (
                                os.getenv("MNPFLEX_PASSWORD") or os.getenv("EPIGNOSTIX_PASSWORD")
                            )
                            if _mnpflex_username and _mnpflex_password:
                                self.bulk_mnpflex_button = ui.button(
                                    "mnpflex run all",
                                    on_click=lambda: None,
                                ).props(
                                    "color=secondary dense no-caps"
                                ).classes(
                                    "border border-slate-300 dark:border-slate-600"
                                )
                            else:
                                self.bulk_mnpflex_button = None

                            self.export_reports_button = ui.button(
                                "Export reports",
                                on_click=lambda: None,
                            ).props("color=primary").classes(
                                "rounded-lg px-4 text-title-medium"
                            )
                            self.export_reports_button.disable()
                            logging.info(
                                "[samples_overview] toolbar controls created "
                                "(Select all / Clear / SNP / Export)"
                            )

                    # Filters row (search left, origin right — read .value in handlers; Quasar
                    # update:model-value payloads are not reliably parsed by _extract_event_value)
                    with ui.row().classes(
                        "w-full items-center justify-between gap-2 mb-2 flex-wrap "
                        "text-slate-700 dark:text-slate-300"
                    ):

                        # Initialize filters model
                        self._samples_filters = getattr(
                            self, "_samples_filters", None
                        ) or {"query": "", "origin": "All", "job_type": "All"}

                        # Global search box
                        self.samples_search = (
                            ui.input(placeholder="Search…")
                            .props("clearable dense outlined")
                            .on("change", self._on_samples_search_update)
                            .on("update:model-value", self._on_samples_search_update)
                            .classes("min-w-[12rem] flex-1 max-w-xl")
                        )

                        # Initialize search input with current filter value
                        try:
                            current_query = self._samples_filters.get("query", "")
                            if current_query:
                                self.samples_search.value = current_query
                        except Exception as e:
                            logging.debug(
                                "[samples_overview] could not preset search filter: %s", e
                            )

                        # Origin filter
                        self.origin_filter = (
                            ui.select(
                                options=["All", "Live", "Pre-existing", "Complete"],
                                value=self._samples_filters.get("origin", "All"),
                                label="Origin",
                            )
                            # Do not set Quasar emit-value — it changes event shape and breaks
                            # NiceGUI Select._event_args_to_value (expects dict args, not raw index).
                            .props("dense clearable outlined")
                            .on("change", self._on_origin_filter_update)
                            .on("update:model-value", self._on_origin_filter_update)
                            .classes("min-w-[11rem]")
                        )

                        # Job type filter (options populated from current row data)
                        self.job_type_filter = (
                            ui.select(
                                options=["All"],
                                value=self._samples_filters.get("job_type", "All"),
                                label="Job type",
                            )
                            .props("dense clearable outlined")
                            .on("change", self._on_job_type_filter_update)
                            .on("update:model-value", self._on_job_type_filter_update)
                            .classes("min-w-[12rem]")
                        )

                    # Loading state container
                    self.samples_loading_container = ui.card().classes(
                        "w-full p-6 text-center border border-slate-200 dark:border-slate-800 "
                        "bg-gradient-to-b from-slate-50 to-white dark:from-slate-900/80 dark:to-zinc-950/90 "
                        "rounded-xl"
                    )
                    with self.samples_loading_container:
                        ui.spinner(size="lg", color="primary")
                        ui.label("Loading samples…").classes(
                            "ml-2 text-title-medium text-slate-800 dark:text-slate-100"
                        )
                        ui.label("This may take a moment for large directories").classes(
                            "text-body-small text-slate-500 dark:text-slate-400 mt-2"
                        )
                    
                    # Create a placeholder table that will be updated later
                    from robin.gui.theme import styled_table

                    # Create samples table
                    _samples_container, self.samples_table = styled_table(
                        columns=[
                            {
                                "name": "actions",
                                "label": "Actions",
                                "field": "actions",
                            },
                            {
                                "name": "sample_id",
                                "label": "Library ID",
                                "field": "sample_id",
                                "sortable": True,
                            },
                            {
                                "name": "test_id",
                                "label": "Test ID",
                                "field": "test_id",
                                "sortable": True,
                            },
                            {
                                "name": "origin",
                                "label": "Origin",
                                "field": "origin",
                                "sortable": True,
                            },
                            {
                                "name": "run_start",
                                "label": "Run Start",
                                "field": "run_start",
                                "sortable": True,
                            },
                            {
                                "name": "device",
                                "label": "Device",
                                "field": "device",
                                "sortable": True,
                            },
                            {
                                "name": "flowcell",
                                "label": "Flowcell",
                                "field": "flowcell",
                                "sortable": True,
                            },
                            {
                                "name": "file_progress",
                                "label": "Job Progress",
                                "field": "file_progress",
                                "sortable": True,
                            },
                            {
                                "name": "pipeline_progress",
                                "label": "Finalize/SNP",
                                "field": "pipeline_progress",
                                "sortable": True,
                            },
                            {
                                "name": "active_jobs",
                                "label": "A",
                                "field": "active_jobs",
                                "sortable": True,
                                "align": "center",
                                "style": "width:52px; max-width:52px;",
                                "headerStyle": "width:52px; max-width:52px;",
                            },
                            {
                                "name": "pending_jobs",
                                "label": "P",
                                "field": "pending_jobs",
                                "sortable": True,
                                "align": "center",
                                "style": "width:52px; max-width:52px;",
                                "headerStyle": "width:52px; max-width:52px;",
                            },
                            {
                                "name": "total_jobs",
                                "label": "T",
                                "field": "total_jobs",
                                "sortable": True,
                                "align": "center",
                                "style": "width:52px; max-width:52px;",
                                "headerStyle": "width:52px; max-width:52px;",
                            },
                            {
                                "name": "completed_jobs",
                                "label": "C",
                                "field": "completed_jobs",
                                "sortable": True,
                                "align": "center",
                                "style": "width:52px; max-width:52px;",
                                "headerStyle": "width:52px; max-width:52px;",
                            },
                            {
                                "name": "failed_jobs",
                                "label": "F",
                                "field": "failed_jobs",
                                "sortable": True,
                                "align": "center",
                                "style": "width:52px; max-width:52px;",
                                "headerStyle": "width:52px; max-width:52px;",
                            },
                            {
                                "name": "job_types",
                                "label": "Job Types",
                                "field": "job_types",
                                "sortable": True,
                            },
                            {
                                "name": "last_seen",
                                "label": "Last Activity",
                                "field": "last_seen",
                                "sortable": True,
                            },
                            {
                                "name": "export",
                                "label": "Export",
                                "field": "export",
                            },
                        ],
                        rows=[],
                        pagination=20,
                        class_size="table-xs",
                        row_key="sample_id",
                    )
                    _export_col = any(
                        (c.get("name") == "export") for c in self.samples_table.columns
                    )
                    logging.info(
                        "[samples_overview] styled_table ready: %d columns, export_col=%s",
                        len(self.samples_table.columns),
                        _export_col,
                    )
                    try:
                        self.samples_table.props("rows-per-page-options=[10,20,50,0]")
                    except Exception as e:
                        logging.debug(
                            "[samples_overview] rows-per-page-options props: %s", e
                        )
                    
                    # Set default sorting by last activity in reverse chronological order
                    try:
                        # Try multiple approaches to set default sorting
                        self.samples_table.props("default-sort=last_seen desc")
                        # Alternative approach using sort property
                        self.samples_table.props("sort=last_seen desc")
                        # Set initial sort on the column
                        for col in self.samples_table.columns:
                            if col.get("field") == "last_seen":
                                col["sort"] = "desc"
                                break
                    except Exception as e:
                        logging.debug("[samples_overview] default sort props: %s", e)
                    
                    # Populate table from master record if available (fast initial load)
                    try:
                        if self._samples_master_record:
                            logging.info("Populating table from master record cache")
                            self._refresh_table_from_master()
                    except Exception as e:
                        logging.debug(f"Error populating table from master record: {e}")

                    # Per-row action buttons (Finalize only shows for Complete/Pre-existing samples)
                    try:
                        self.samples_table.add_slot(
                            "body-cell-actions",
                            """
<q-td key=\"actions\" :props=\"props\">
  <q-btn-group>
    <q-btn color=\"primary\" size=\"sm\" label=\"View\"
           :href=\"'/live_data/' + encodeURIComponent(props.row.sample_id)\" />
    <q-btn v-if=\"props.row.origin === 'Complete' || props.row.origin === 'Pre-existing'\"
           color=\"secondary\" size=\"sm\" label=\"Finalize\" icon=\"merge_type\"
           @click=\"$parent.$emit('finalize-target', props.row.sample_id)\" />
  </q-btn-group>
</q-td>
""",
                        )
                    except Exception as e:
                        logging.warning(
                            "[samples_overview] add_slot body-cell-actions failed: %s",
                            e,
                            exc_info=True,
                        )
                    
                    # Add job progress column with linear progress bar
                    try:
                        self.samples_table.add_slot(
                            "body-cell-file_progress",
                            """
<q-td key=\"file_progress\" :props=\"props\">
  <div style=\"min-width: 120px;\">
    <q-linear-progress 
      :value=\"props.row.file_progress || 0\" 
      :color=\"props.row.file_progress >= 1 ? 'positive' : 'primary'\" 
      size=\"12px\" 
      rounded
      class=\"q-mb-xs\" />
    <div class=\"text-center text-[10px] text-slate-600 dark:text-slate-400\">
      {{ props.row.files_processed || 0 }}/{{ props.row.files_seen || 0 }} jobs
    </div>
  </div>
</q-td>
""",
                        )
                    except Exception as e:
                        logging.warning(
                            "[samples_overview] add_slot body-cell-file_progress failed: %s",
                            e,
                            exc_info=True,
                        )

                    # Finalize/SNP progress indicator per sample
                    try:
                        self.samples_table.add_slot(
                            "body-cell-pipeline_progress",
                            """
<q-td key=\"pipeline_progress\" :props=\"props\">
  <!-- Compact: small state dot + one-line phase (ellipsis) -->
  <div class=\"flex items-center\" style=\"gap: 8px; max-width: 170px;\">
    <span
      :style=\"{
        width: '10px',
        height: '10px',
        borderRadius: '999px',
        backgroundColor:
          (props.row.pipeline_status || '').toLowerCase().includes('fail') ? '#f87171' :
          ((props.row.pipeline_progress || 0) >= 1) ? '#22c55e' :
          ((props.row.pipeline_status || '') ? '#60a5fa' : '#94a3b8'),
        boxShadow: '0 0 0 2px rgba(148,163,184,0.15)'
      }\" />
    <div
      class=\"text-[10px] text-slate-700 dark:text-slate-200\"
      style=\"white-space: nowrap; overflow: hidden; text-overflow: ellipsis;\"
      :title=\"props.row.pipeline_detail ? (props.row.pipeline_status + ' - ' + props.row.pipeline_detail) : (props.row.pipeline_status || 'Idle')\"
    >
      {{ props.row.pipeline_status || 'Idle' }}
    </div>
  </div>
</q-td>
""",
                        )
                    except Exception as e:
                        logging.warning(
                            "[samples_overview] add_slot body-cell-pipeline_progress failed: %s",
                            e,
                            exc_info=True,
                        )

                    # Job types: allow wrapping (comma+space separated in row data)
                    try:
                        self.samples_table.add_slot(
                            "body-cell-job_types",
                            """
<q-td key="job_types" :props="props">
  <div class="text-xs whitespace-normal break-words max-w-[min(28rem,50vw)]"
       style="word-break: break-word;">
    {{ props.row.job_types }}
  </div>
</q-td>
""",
                        )
                    except Exception as e:
                        logging.warning(
                            "[samples_overview] add_slot body-cell-job_types failed: %s",
                            e,
                            exc_info=True,
                        )
                    
                    # Add export checkbox as rightmost column + header "select all" checkbox
                    try:
                        self.samples_table.add_slot(
                            "header-cell-export",
                            """
<q-th :props="props">
  <div class="column items-center q-gutter-xs">
    <span class="text-caption text-slate-500">All</span>
    <q-checkbox
      dense
      size="sm"
      @update:model-value="$parent.$emit('export-header-toggle', $event)"
    />
  </div>
</q-th>
""",
                        )
                        self.samples_table.add_slot(
                            "body-cell-export",
                            """
<q-td key=\"export\" :props=\"props\">
  <q-checkbox size=\"sm\"
              :model-value=\"props.row.export === true\"
              @update:model-value=\"$parent.$emit('export-toggled', { id: props.row.sample_id, value: $event })\" />
</q-td>
""",
                        )
                    except Exception as e:
                        logging.warning(
                            "[samples_overview] add_slot export column failed: %s",
                            e,
                            exc_info=True,
                        )
                    
                    # Handle finalize-target action
                    def _on_finalize_target(event):
                        try:
                            logging.debug("finalize-target event: %r", event)
                            sample_id = getattr(event, "args", None) if hasattr(event, "args") else None
                            logging.debug("finalize-target sample_id=%r", sample_id)
                            if isinstance(sample_id, str):
                                # Keep this UI callback non-blocking.
                                # All heavyweight checks/submissions run in background logic.
                                self._trigger_target_bam_finalization(sample_id, trigger_snp=True)
                                logging.debug(
                                    "finalize-target queued: sample_id=%s trigger_snp=True",
                                    sample_id,
                                )
                                ui.notify(
                                    f"Triggered run finalization for {sample_id}. "
                                    "SNP analysis will be queued if prerequisites are met.",
                                    type="positive",
                                )
                            else:
                                ui.notify("Invalid sample ID", type="warning")
                        except Exception as e:
                            ui.notify(f"Error triggering finalization: {e}", type="negative")
                    
                    try:
                        self.samples_table.on("finalize-target", _on_finalize_target)
                    except Exception as e:
                        logging.warning(
                            "[samples_overview] bind finalize-target failed: %s",
                            e,
                            exc_info=True,
                        )

                    # Track multi-selection for batch export via custom checkbox column
                    try:
                        self._selected_sample_ids = set()

                        def _on_export_toggled(event):
                            try:
                                payload = None
                                if hasattr(event, "args"):
                                    payload = getattr(event, "args", None)
                                elif isinstance(event, dict):
                                    payload = event
                                if isinstance(payload, dict):
                                    sid = payload.get("id")
                                    val = bool(payload.get("value"))
                                    if sid:
                                        if val:
                                            self._selected_sample_ids.add(sid)
                                        else:
                                            self._selected_sample_ids.discard(sid)
                                        # reflect state back into rows
                                        try:
                                            for r in self.samples_table.rows or []:
                                                if r.get("sample_id") == sid:
                                                    r["export"] = (
                                                        sid in self._selected_sample_ids
                                                    )
                                            self.samples_table.update()
                                        except Exception as ue:
                                            logging.debug(
                                                "[samples_overview] export row sync failed: %s",
                                                ue,
                                                exc_info=True,
                                            )
                                        if self._selected_sample_ids:
                                            self.export_reports_button.enable()
                                        else:
                                            self.export_reports_button.disable()
                                        logging.info(
                                            "[samples_overview] export-toggled id=%s value=%s "
                                            "selected_count=%s",
                                            sid,
                                            val,
                                            len(self._selected_sample_ids),
                                        )
                            except Exception as e:
                                logging.warning(
                                    "[samples_overview] export-toggled handler error: %s",
                                    e,
                                    exc_info=True,
                                )

                        self.samples_table.on("export-toggled", _on_export_toggled)

                        def _on_export_header_toggle(event):
                            try:
                                logging.info(
                                    "[samples_overview] export-header-toggle raw_args=%r",
                                    getattr(event, "args", None)
                                    if hasattr(event, "args")
                                    else event,
                                )
                                val = True
                                if hasattr(event, "args"):
                                    a = getattr(event, "args", None)
                                    if isinstance(a, (list, tuple)) and len(a) > 0:
                                        val = bool(a[0])
                                    elif isinstance(a, dict):
                                        val = bool(a.get("value", True))
                                    else:
                                        val = bool(a)
                                elif isinstance(event, dict):
                                    val = bool(event.get("value", True))
                                if val:
                                    self._samples_select_all_visible_for_export()
                                else:
                                    self._samples_clear_export_selection()
                            except Exception as e:
                                logging.warning(
                                    "[samples_overview] export-header-toggle failed (%s); "
                                    "falling back to select-all",
                                    e,
                                    exc_info=True,
                                )
                                self._samples_select_all_visible_for_export()

                        self.samples_table.on(
                            "export-header-toggle", _on_export_header_toggle
                        )
                        logging.info(
                            "[samples_overview] export-toggled and export-header-toggle "
                            "handlers registered"
                        )
                    except Exception as e:
                        logging.warning(
                            "[samples_overview] export checkbox table.on wiring failed: %s",
                            e,
                            exc_info=True,
                        )

                    # Batch export selected reports + wire SNP bulk (separate try so one failure does not block the other)
                    try:
                        logging.info(
                            "[samples_overview] wiring bulk_snp and export_reports "
                            "on_click handlers"
                        )

                        async def _ask_run_bulk_snp():
                            logging.info("[samples_overview] SNP: all missing clicked")
                            logging.info(
                                "[samples_overview] SNP click handler thread=%s",
                                threading.current_thread().name,
                            )
                            if self._bulk_snp_worker_running:
                                ui.notify(
                                    "SNP batch is already running.",
                                    type="warning",
                                )
                                return
                            t0 = time.time()
                            # Offload filesystem scan so we don't block the NiceGUI/UI event loop.
                            missing = await asyncio.to_thread(
                                self._list_samples_needing_snp_calling
                            )
                            logging.info(
                                "[samples_overview] SNP missing scan done in %.2fs (n=%d)",
                                time.time() - t0,
                                len(missing or []),
                            )
                            if not missing:
                                ui.notify(
                                    "No samples need SNP calling (outputs exist or prerequisites missing).",
                                    type="info",
                                )
                                return
                            n = len(missing)
                            with ui.dialog() as dlg:
                                with ui.card().classes(
                                    "robin-dialog-surface w-96 max-w-[95vw] p-4"
                                ):
                                    ui.label("Run SNP for all missing samples?").classes(
                                        "classification-insight-heading text-headline-small mb-2"
                                    )
                                    ui.label(
                                        f"{n} sample(s) will be processed one after another "
                                        "(same as triggering SNP manually each time). "
                                        "This may take a long time overall."
                                    ).classes("text-sm text-gray-600 mb-4")
                                    with ui.row().classes("justify-end gap-2 flex-wrap"):
                                        ui.button(
                                            "Cancel",
                                            on_click=lambda: dlg.submit(False),
                                        ).props("flat no-caps outline")
                                        ui.button(
                                            "Start",
                                            on_click=lambda: dlg.submit(True),
                                        ).props("color=primary no-caps")
                            if not await dlg:
                                return
                            self._bulk_snp_worker_running = True
                            ids_copy = list(missing)

                            def _run():
                                self._bulk_snp_sequential_run(ids_copy)

                            threading.Thread(
                                target=_run,
                                daemon=True,
                                name="robin-bulk-snp",
                            ).start()
                            ui.notify(
                                f"Started sequential SNP batch for {n} sample(s).",
                                type="positive",
                            )

                        try:
                            self.bulk_snp_button.on_click(_ask_run_bulk_snp)
                        except Exception as e:
                            logging.warning(
                                "Could not wire SNP bulk button (export still works): %s",
                                e,
                            )

                        # Optional MNP-Flex bulk run (same pattern as SNP)
                        try:
                            if getattr(self, "bulk_mnpflex_button", None) is not None:

                                async def _ask_run_bulk_mnpflex():
                                    logging.info(
                                        "[samples_overview] mnpflex run all clicked"
                                    )
                                    logging.info(
                                        "[samples_overview] mnpflex click handler thread=%s",
                                        threading.current_thread().name,
                                    )
                                    if self._bulk_mnpflex_worker_running:
                                        ui.notify(
                                            "MNP-Flex batch is already running.",
                                            type="warning",
                                        )
                                        return
                                    if not self._is_mnpflex_enabled_for_gui():
                                        ui.notify(
                                            "MNP-Flex is not enabled (missing credentials).",
                                            type="warning",
                                        )
                                        return

                                    t0 = time.time()
                                    # Offload disk checks so the UI thread stays responsive.
                                    rows_override = list(
                                        getattr(self.samples_table, "rows", None)
                                        or getattr(self, "_last_samples_rows", [])
                                        or []
                                    )
                                    missing = await asyncio.to_thread(
                                        self._list_samples_needing_mnpflex_for_visible,
                                        rows_override,
                                    )
                                    logging.info(
                                        "[samples_overview] mnpflex missing scan done in %.2fs (n=%d)",
                                        time.time() - t0,
                                        len(missing or []),
                                    )
                                    if not missing:
                                        ui.notify(
                                            "No samples need MNP-Flex (results exist or prerequisites missing).",
                                            type="info",
                                        )
                                        return
                                    n = len(missing)

                                    with ui.dialog() as dlg:
                                        with ui.card().classes(
                                            "robin-dialog-surface w-96 max-w-[95vw] p-4"
                                        ):
                                            ui.label(
                                                "Run MNP-Flex for eligible missing samples?"
                                            ).classes(
                                                "classification-insight-heading text-headline-small mb-2"
                                            )
                                            ui.label(
                                                f"{n} sample(s) will be processed one after another. "
                                                "This can take a long time overall."
                                            ).classes(
                                                "text-sm text-gray-600 mb-4"
                                            )
                                            with ui.row().classes(
                                                "justify-end gap-2 flex-wrap"
                                            ):
                                                ui.button(
                                                    "Cancel",
                                                    on_click=lambda: dlg.submit(False),
                                                ).props(
                                                    "flat no-caps outline"
                                                )
                                                ui.button(
                                                    "Start",
                                                    on_click=lambda: dlg.submit(True),
                                                ).props("color=primary no-caps")

                                    if not await dlg:
                                        return

                                    self._bulk_mnpflex_worker_running = True
                                    ids_copy = list(missing)

                                    def _run():
                                        self._bulk_mnpflex_sequential_run(ids_copy)

                                    threading.Thread(
                                        target=_run,
                                        daemon=True,
                                        name="robin-bulk-mnpflex",
                                    ).start()
                                    ui.notify(
                                        f"Started sequential MNP-Flex batch for {n} sample(s).",
                                        type="positive",
                                    )

                                self.bulk_mnpflex_button.on_click(_ask_run_bulk_mnpflex)
                        except Exception as e:
                            logging.warning(
                                "Could not wire MNP-Flex bulk button: %s",
                                e,
                                exc_info=True,
                            )

                        async def _export_selected_reports(state: Dict[str, Any], progress_dialog, files_to_download, download_complete, progress_callback, progress_updates):
                            try:
                                selected = list(
                                    getattr(self, "_selected_sample_ids", set()) or []
                                )
                                if not selected:
                                    ui.notify("No samples selected", type="warning")
                                    return
                                
                                total_samples = len(selected)
                                
                                try:
                                    from nicegui import run as ng_run  # type: ignore
                                except Exception:
                                    ng_run = None  # type: ignore
                                
                                for idx, sid in enumerate(selected):
                                    try:
                                        # Update overall progress
                                        overall_progress = (idx / total_samples) * 0.9
                                        
                                        # Emit progress update showing current sample and mark as starting
                                        current_sample_msg = f"Generating {idx + 1}/{total_samples} - {sid}"
                                        progress_updates.put({
                                            'stage': 'processing_sections',
                                            'message': 'Starting...',
                                            'progress': 0.0,
                                            'sample_id': sid
                                        })
                                        
                                        sample_dir = (
                                            Path(self.monitored_directory) / sid
                                            if self.monitored_directory
                                            else None
                                        )
                                        if not sample_dir or not sample_dir.exists():
                                            logging.warning(f"Missing output for {sid}")
                                            continue
                                        
                                        # Don't use notification system - only update dialog
                                        
                                        filename = f"{sid}_run_report.pdf"
                                        pdf_path = os.path.join(
                                            str(sample_dir), filename
                                        )
                                        os.makedirs(str(sample_dir), exist_ok=True)
                                        export_csv_dir = None
                                        if bool(state.get("export_csv", False)):
                                            export_csv_dir = os.path.join(
                                                str(sample_dir), "report_csv"
                                            )
                                        
                                        # Don't use the notification system - use only our dialog callback
                                        if ng_run is not None:
                                            # Use custom callback that updates dialog only
                                            def sample_progress_callback(data: Dict[str, Any]):
                                                data['sample_id'] = sid  # Add sample ID to track which bar to update
                                                progress_callback(data)  # This updates the dialog via progress_updates queue
                                            
                                            pdf_file = await ng_run.io_bound(
                                                create_pdf,
                                                pdf_path,
                                                str(sample_dir),
                                                self.center or "Unknown",
                                                report_type=state.get("type", "detailed"),
                                                export_csv_dir=export_csv_dir,
                                                export_xlsx=False,
                                                export_zip=bool(
                                                    state.get("export_csv", False)
                                                ),
                                                progress_callback=sample_progress_callback,  # Pass the dialog-only callback
                                                workflow_steps=self.workflow_steps if hasattr(self, 'workflow_steps') else None,
                                            )
                                        else:
                                            # Use custom callback that updates dialog only
                                            def sample_progress_callback(data: Dict[str, Any]):
                                                data['sample_id'] = sid  # Add sample ID to track which bar to update
                                                progress_callback(data)  # This updates the dialog via progress_updates queue
                                            
                                            pdf_file = create_pdf(
                                                pdf_path,
                                                str(sample_dir),
                                                self.center or "Unknown",
                                                report_type=state.get("type", "detailed"),
                                                export_csv_dir=export_csv_dir,
                                                export_xlsx=False,
                                                export_zip=bool(
                                                    state.get("export_csv", False)
                                                ),
                                                progress_callback=sample_progress_callback,  # Pass the dialog-only callback
                                                workflow_steps=self.workflow_steps if hasattr(self, 'workflow_steps') else None,
                                            )
                                        
                                        # Queue file for download instead of downloading immediately
                                        files_to_download.append(pdf_file)
                                        
                                        # Also offer CSV ZIP if requested
                                        if (
                                            bool(state.get("export_csv", False))
                                            and export_csv_dir
                                        ):
                                            zip_path = os.path.join(
                                                export_csv_dir, f"{sid}_report_data.zip"
                                            )
                                            if os.path.exists(zip_path):
                                                files_to_download.append(zip_path)
                                        
                                        # Mark sample as complete
                                        progress_updates.put({
                                            'stage': 'completed',
                                            'message': 'Completed',
                                            'progress': 1.0,
                                            'sample_id': sid
                                        })
                                        
                                    except Exception as e:
                                        # Report generation failed
                                        logging.error(f"Export failed for {sid}: {e}")
                                        # Mark sample as failed
                                        progress_updates.put({
                                            'stage': 'error',
                                            'message': f'Failed: {str(e)[:50]}',
                                            'progress': 1.0,
                                            'sample_id': sid
                                        })
                                
                                # Mark as complete
                                download_complete["done"] = True
                                logging.info(f"Bulk export complete. {len(files_to_download)} file(s) ready for download.")
                            except Exception as e:
                                logging.error(f"Error in bulk export: {e}")
                                download_complete["done"] = True

                        async def _confirm_bulk_export():
                            logging.info(
                                "[samples_overview] Export reports clicked "
                                "selected=%s",
                                len(
                                    getattr(self, "_selected_sample_ids", None) or []
                                ),
                            )
                            # Ensure there is at least one selection before opening
                            if not getattr(self, "_selected_sample_ids", None):
                                ui.notify("No samples selected", type="warning")
                                return

                            selected_ids = list(getattr(self, "_selected_sample_ids", set()) or [])
                            num_selected = len(selected_ids)

                            report_types = {
                                "summary": "Summary Only",
                                "detailed": "Detailed",
                            }
                            state: Dict[str, Any] = {
                                "type": "detailed",
                                "export_csv": False,
                            }

                            with ui.dialog().props("persistent") as dialog:
                                with ui.card().classes(
                                    "robin-dialog-surface w-96 max-w-[95vw] p-4"
                                ):
                                    ui.label("Export reports").classes(
                                        "classification-insight-heading text-headline-small mb-3"
                                    )

                                    with ui.column():
                                        with ui.column().classes("mb-4"):
                                            ui.label("Report type").classes(
                                                "target-coverage-panel__meta-label mb-2"
                                            )
                                            ui.toggle(
                                                report_types,
                                                value="detailed",
                                                on_change=lambda e: state.update(
                                                    {"type": e.value}
                                                ),
                                            )

                                        with ui.column().classes("mb-4"):
                                            ui.label("Include data").classes(
                                                "target-coverage-panel__meta-label mb-2"
                                            )
                                            ui.checkbox(
                                                "CSV data (ZIP)",
                                                value=False,
                                                on_change=lambda e: state.update(
                                                    {"export_csv": bool(e.value)}
                                                ),
                                            )

                                        with ui.column().classes("mb-4"):
                                            ui.label("Disclaimer").classes(
                                                "target-coverage-panel__meta-label mb-2"
                                            )
                                            formatted_text = EXTENDED_DISCLAIMER_TEXT.replace(
                                                "\n\n", "<br><br>"
                                            ).replace("\n", " ")
                                            ui.label(formatted_text).classes(
                                                "text-sm text-gray-600 mb-4"
                                            )

                                        ui.label(
                                            f"Are you sure you want to export reports for {num_selected} sample(s)?"
                                        ).classes("classification-insight-foot mb-4")

                                        with ui.row().classes("justify-end gap-2 flex-wrap"):
                                            ui.button(
                                                "Cancel",
                                                on_click=lambda: dialog.submit("Cancel"),
                                            ).props("flat no-caps outline")
                                            ui.button(
                                                "Export",
                                                on_click=lambda: dialog.submit("Export"),
                                                icon="download",
                                            ).props("color=primary no-caps")

                                    ui.label(
                                        f"Exporting reports for {num_selected} sample(s): {', '.join(selected_ids[:3])}"
                                        + (
                                            f" and {num_selected - 3} more"
                                            if num_selected > 3
                                            else ""
                                        )
                                    ).classes(
                                        "text-sm font-medium text-gray-700 mt-4"
                                    )

                            dialog_result = await dialog
                            if dialog_result != "Export":
                                return

                            # Now show the progress dialog
                            with ui.dialog().props("persistent") as progress_dialog:
                                with ui.card().classes(
                                    "robin-dialog-surface w-96 max-w-[95vw] p-4"
                                ):
                                    ui.label("Exporting reports").classes(
                                        "classification-insight-heading text-headline-small mb-3"
                                    )

                                    ui.label(f"Exporting {num_selected} report(s)").classes(
                                        "text-sm font-medium text-gray-700 mb-4"
                                    )

                                    # Report type and output displays
                                    ui.label(
                                        f"Report Type: {state.get('type', 'detailed').title()}"
                                    ).classes("text-sm mb-2")
                                    
                                    ui.label(
                                        f"Output: PDF{' + CSV (ZIP)' if state.get('export_csv', False) else ''}"
                                    ).classes("text-sm mb-4")
                                    
                                    # Create individual progress bars for each sample
                                    sample_progress_bars = {}
                                    sample_progress_labels = {}
                                    
                                    with ui.column().classes("w-full"):
                                        for sid in selected_ids:
                                            with ui.column().classes("mb-3 w-full"):
                                                ui.label(sid).classes("text-xs font-medium text-gray-700 mb-1")
                                                progress_bar = ui.linear_progress(0.0).classes("mb-1")
                                                progress_label = ui.label("Waiting...").classes("text-xs text-gray-500")
                                                sample_progress_bars[sid] = progress_bar
                                                sample_progress_labels[sid] = progress_label
                                    
                                    # Messages container - use label with newlines for multiple messages
                                    messages_label = ui.label("").classes("text-xs text-gray-500")

                                    # Track messages
                                    progress_updates = queue.Queue()
                                    messages_list = []
                                    
                                    # Track current sample being processed
                                    current_sample = {"id": None}
                                    
                                    # Timer to process progress updates on UI thread
                                    def process_progress_updates():
                                        """Process queued progress updates."""
                                        try:
                                            while not progress_updates.empty():
                                                update = progress_updates.get_nowait()
                                                stage = update.get("stage", "unknown")
                                                message = update.get("message", "")
                                                progress = update.get("progress", 0.0)
                                                sample_id = update.get("sample_id")
                                                
                                                # Update the current sample's progress bar
                                                if sample_id and sample_id in sample_progress_bars:
                                                    if progress is not None:
                                                        sample_progress_bars[sample_id].value = progress
                                                        sample_progress_labels[sample_id].text = f"{int(progress * 100)}% - {message}"
                                                    else:
                                                        sample_progress_labels[sample_id].text = message
                                                    current_sample["id"] = sample_id
                                                
                                                # Add message to messages list
                                                messages_list.append(message)
                                                # Keep only last 10 messages
                                                if len(messages_list) > 10:
                                                    messages_list.pop(0)
                                                
                                                # Update messages label
                                                messages_label.text = "\n".join(messages_list[-5:])
                                                
                                        except queue.Empty:
                                            pass
                                        except Exception as e:
                                            logging.debug(f"Error processing progress updates: {e}")
                                    
                                    # Set up timer to process updates
                                    update_timer = ui.timer(0.1, process_progress_updates)
                                    
                                    def progress_callback(progress_data: Dict[str, Any]):
                                        """Custom progress callback to update dialog (called from background thread)."""
                                        try:
                                            # Queue the update instead of directly updating UI
                                            progress_updates.put(progress_data)
                                        except Exception as e:
                                            logging.debug(f"Error in progress callback: {e}")
                                    
                                    # Track if still generating
                                    is_generating = {"active": True}

                                    # Storage for the files to download
                                    files_to_download = []
                                    download_complete = {"done": False}
                                    
                                    # Timer to handle downloads once background task is done
                                    def handle_downloads():
                                        """Handle downloads in UI context once generation is complete."""
                                        if download_complete["done"] and files_to_download:
                                            # Log that export is complete
                                            valid_files = [
                                                f
                                                for f in files_to_download
                                                if f is not None and os.path.isfile(f)
                                            ]
                                            logging.info(
                                                f"Bulk export complete. {len(valid_files)} file(s) ready for download."
                                            )
                                            # Safari blocks multiple programmatic downloads; use one ZIP when needed.
                                            bundle = self._zip_paths_for_bulk_download(valid_files)
                                            if bundle:
                                                if len(valid_files) > 1:
                                                    ui.notify(
                                                        "Downloading a single ZIP (works in Safari; multiple separate downloads are blocked there).",
                                                        type="info",
                                                    )
                                                logging.debug(f"Initiating download: {bundle}")
                                                ui.download(bundle)
                                                if bundle.endswith(".zip") and "robin_reports_" in os.path.basename(
                                                    bundle
                                                ):
                                                    ui.timer(
                                                        120.0,
                                                        lambda p=bundle: self._unlink_quiet(p),
                                                        once=True,
                                                    )
                                            # Close dialog after 3 seconds
                                            ui.timer(3.0, lambda: progress_dialog.submit(None), once=True)
                                            download_timer.deactivate()
                                    
                                    download_timer = ui.timer(0.1, handle_downloads)
                                    
                                    # Start report generation
                                    async def complete_export():
                                        """Complete the export and close dialog."""
                                        try:
                                            await _export_selected_reports(state, progress_dialog, files_to_download, download_complete, progress_callback, progress_updates)
                                        finally:
                                            # Clean up timer
                                            update_timer.deactivate()
                                            is_generating["active"] = False
                                    
                                    asyncio.create_task(complete_export())

                            await progress_dialog

                        # Wire the button now that handlers exist
                        self.export_reports_button.on_click(_confirm_bulk_export)
                        logging.info(
                            "[samples_overview] bulk SNP and Export reports buttons "
                            "on_click wired successfully"
                        )
                    except Exception:
                        logging.exception(
                            "[samples_overview] bulk export / SNP wiring failed — "
                            "Export and SNP buttons may still use placeholder handlers"
                        )

                    # If we have cached rows, show them immediately for instant UX
                    if self._last_samples_rows:
                        try:
                            self._apply_samples_table_filters()
                            # Hide loading container if we have data
                            if hasattr(self, "samples_loading_container"):
                                self.samples_loading_container.set_visibility(False)
                            # No selection behavior needed; per-row buttons handle navigation
                        except Exception as e:
                            logging.debug(
                                "[samples_overview] apply cached rows on load: %s", e
                            )
                    else:
                        # Show loading container if no cached data
                        if hasattr(self, "samples_loading_container"):
                            self.samples_loading_container.set_visibility(True)


    def _update_samples_table_sync(self, data: Dict[str, Any]) -> None:
        """Synchronous version of samples table update"""
        try:
            if not hasattr(self, "samples_table"):
                return
            samples = data.get("samples", [])
            expected_job_types = self._get_expected_completion_job_types()
            base = (
                Path(self.monitored_directory) if self.monitored_directory else None
            )

            # Deduplicate by sample_id taking the newest last_seen
            by_id: Dict[str, Dict[str, Any]] = {}
            for s in samples:
                sid = s.get("sample_id", "") or "unknown"
                last_seen = float(s.get("last_seen", time.time()))
                existing = by_id.get(sid)
                if not existing or last_seen >= existing.get("_last_seen_raw", 0):
                    origin_value = (
                        "Pre-existing"
                        if sid in self._preexisting_sample_ids and (time.time() - last_seen) >= self.completion_timeout_seconds
                        else "Live"
                    )
                    # Flip Live samples to Complete if inactive for configured timeout
                    prev_origin = None
                    if existing:
                        prev_origin = existing.get("origin")
                    else:
                        # Check previous rows for this sample
                        existing_rows = self._last_samples_rows or []
                        for prev_row in existing_rows:
                            if prev_row.get("sample_id") == sid:
                                prev_origin = prev_row.get("origin")
                                break
                    merged = self._merge_workflow_sample_job_counts(sid, s)
                    if (
                        origin_value == "Live"
                        and base is not None
                        and (time.time() - last_seen) >= self.completion_timeout_seconds
                    ):
                        if merged["active_jobs"] == 0 and merged["pending_jobs"] == 0:
                            sample_dir = base / sid
                            if expected_job_types:
                                complete_on_disk = self._expected_jobs_completed(
                                    sample_dir, expected_job_types
                                )
                            else:
                                complete_on_disk = self._is_target_bam_finalize_redundant(
                                    sid
                                )
                            if not complete_on_disk:
                                complete_on_disk = self._is_target_bam_finalize_redundant(
                                    sid
                                )
                            if complete_on_disk:
                                origin_value = "Complete"
                    try:
                        active_jobs_count = merged["active_jobs"]
                        pending_jobs_count = merged["pending_jobs"]
                        # Only mark as Complete if timeout passed AND no active jobs
                        if origin_value == "Live" and (time.time() - last_seen) >= self.completion_timeout_seconds:
                            if active_jobs_count == 0 and pending_jobs_count == 0:
                                should_complete = True
                                if expected_job_types and base is not None:
                                    sample_dir = base / sid
                                    should_complete = self._expected_jobs_completed(
                                        sample_dir, expected_job_types
                                    )
                                if should_complete:
                                    origin_value = "Complete"
                                    # Trigger finalization only on a real Live→Complete transition
                                    if prev_origin == "Live":
                                        self._trigger_target_bam_finalization(sid)
                            # If there are active jobs, keep as Live even if timeout passed
                    except Exception:
                        pass
                    total_jobs = merged["total_jobs"]
                    completed_jobs = merged["completed_jobs"]
                    failed_jobs = merged["failed_jobs"]
                    
                    # Set file progress directly from job counts (same data source as other columns)
                    files_seen = total_jobs
                    files_processed = completed_jobs
                    file_progress = completed_jobs / total_jobs if total_jobs > 0 else 0.0
                    
                    by_id[sid] = {
                        "sample_id": sid,
                        "origin": origin_value,
                        # Persisted run info will be patched in below from master.csv
                        "run_start": "",
                        "device": "",
                        "flowcell": "",
                        "active_jobs": merged["active_jobs"],
                        "pending_jobs": merged["pending_jobs"],
                        "total_jobs": total_jobs,
                        "completed_jobs": completed_jobs,
                        "failed_jobs": failed_jobs,
                        "job_types": self._merge_job_types_with_persisted(sid, s),
                        "last_seen": time.strftime(
                            "%Y-%m-%d %H:%M:%S", time.localtime(last_seen)
                        ),
                        "files_seen": files_seen,
                        "files_processed": files_processed,
                        "file_progress": file_progress,
                        "actions": "View",
                        "_last_seen_raw": last_seen,
                    }

            # Patch from master.csv and persist the new overview values for later reload
            try:
                manager = (
                    MasterCSVManager(str(base)) if base and base.exists() else None
                )
            except Exception:
                base, manager = None, None

            for sid, row in by_id.items():
                try:
                    # Persist overview numbers to master.csv so we can restore later
                    if manager is not None:
                        persist_payload = {
                            "active_jobs": int(row.get("active_jobs", 0)),
                            "pending_jobs": int(row.get("pending_jobs", 0)),
                            "total_jobs": int(row.get("total_jobs", 0)),
                            "completed_jobs": int(row.get("completed_jobs", 0)),
                            "failed_jobs": int(row.get("failed_jobs", 0)),
                            "job_types": row.get("job_types", ""),
                            "last_seen": float(row.get("_last_seen_raw", time.time())),
                        }
                        manager.update_sample_overview(sid, persist_payload)

                    # Read run info and file progress from master.csv to display in table
                    if base is not None:
                        csv_path = base / sid / "master.csv"
                        if csv_path.exists():
                            with csv_path.open("r", newline="") as fh:
                                reader = csv.DictReader(fh)
                                first_row = next(reader, None)
                            if first_row:
                                row["run_start"] = self._format_timestamp_for_display(
                                    first_row.get("run_info_run_time", "")
                                )
                                row["device"] = first_row.get("run_info_device", "")
                                row["flowcell"] = first_row.get(
                                    "run_info_flow_cell", ""
                                )
                except Exception:
                    pass

            # Merge with preexisting scans if any
            existing_rows_by_id = {
                r["sample_id"]: r for r in (self._last_samples_rows or [])
            }
            for sid, row in by_id.items():
                existing_rows_by_id[sid] = row
            rows = list(existing_rows_by_id.values())
            # Replace rows to avoid duplicates then apply filters
            self._last_samples_rows = rows
            
            # Update UI on main thread
            if rows and hasattr(self, "samples_table"):
                self._apply_samples_table_filters()
                # Also update file progress display in workflow monitor
                self._update_file_progress()

        except Exception as e:
            logging.error(f"Error in samples table update: {e}")

    # -------- Samples table helpers: search, filter, sort --------
    def _format_timestamp_for_display(self, value: Any) -> str:
        """Convert timestamps (ISO 8601 or epoch) to 'YYYY-MM-DD HH:MM'."""
        try:
            if value is None:
                return ""
            if isinstance(value, (int, float)):
                return datetime.fromtimestamp(float(value)).strftime("%Y-%m-%d %H:%M")
            s = str(value).strip()
            if not s:
                return ""
            if s.endswith("Z"):
                s = s[:-1] + "+00:00"
            try:
                return datetime.fromisoformat(s).strftime("%Y-%m-%d %H:%M")
            except Exception:
                try:
                    return datetime.fromtimestamp(float(s)).strftime("%Y-%m-%d %H:%M")
                except Exception:
                    return s
        except Exception:
            try:
                return str(value)
            except Exception:
                return ""

    def _format_job_types_display(self, job_types: Any) -> str:
        """Normalize job types to comma + space separated text so cells wrap in the table."""
        if job_types is None:
            return ""
        if isinstance(job_types, (list, set, tuple)):
            parts = sorted({str(x).strip() for x in job_types if str(x).strip()})
            return ", ".join(parts)
        s = str(job_types).strip()
        if not s:
            return ""
        parts = [p.strip() for p in s.replace(", ", ",").split(",") if p.strip()]
        return ", ".join(parts)

    def _get_samples_search_query(self) -> str:
        try:
            return (self.samples_search.value or "").strip()
        except Exception:
            return ""

    def _on_samples_search_update(self, _=None) -> None:
        self._set_samples_query(self._get_samples_search_query())

    def _get_origin_filter_value(self) -> str:
        try:
            v = self.origin_filter.value
            if v is None or v == "":
                return "All"
            return str(v)
        except Exception:
            return "All"

    def _on_origin_filter_update(self, _=None) -> None:
        self._set_samples_origin_filter(self._get_origin_filter_value())

    def _get_job_type_filter_value(self) -> str:
        try:
            v = self.job_type_filter.value
            if v is None or v == "":
                return "All"
            return str(v)
        except Exception:
            return "All"

    def _on_job_type_filter_update(self, _=None) -> None:
        self._set_samples_job_type_filter(self._get_job_type_filter_value())

    def _normalize_rows_for_display(
        self, rows: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        normalized: List[Dict[str, Any]] = []
        for r in rows or []:
            sid = r.get("sample_id")
            if not sid or sid == "unknown":
                continue
            jt = r.get("job_types")
            r = dict(r)
            r["job_types"] = self._format_job_types_display(jt)
            normalized.append(r)
        return normalized

    def _extract_event_value(self, event, default=""):
        """Extract value from NiceGUI event arguments."""
        try:
            if event is None:
                return str(default)
            if isinstance(event, (str, int, float, bool)):
                return str(event)
            if hasattr(event, "value"):
                v = event.value
                if v is None:
                    return str(default)
                if isinstance(v, dict):
                    return str(
                        v.get("value", v.get("label", v.get("modelValue", default)))
                    )
                return str(v)
            if hasattr(event, "args"):
                args = event.args
                if args is None:
                    return str(default)
                if isinstance(args, dict):
                    return str(
                        args.get("value", args.get("label", args.get("modelValue", default)))
                    )
                if isinstance(args, (list, tuple)) and len(args) > 0:
                    first = args[0]
                    if isinstance(first, dict):
                        return str(
                            first.get(
                                "value",
                                first.get("label", first.get("modelValue", default)),
                            )
                        )
                    return str(first)
                return str(args)
            return str(default)
        except Exception:
            return str(default)

    def _set_samples_query(self, query: str) -> None:
        try:
            # Ensure we have a filters dict
            if not hasattr(self, "_samples_filters"):
                self._samples_filters = {
                    "query": "",
                    "origin": "All",
                    "job_type": "All",
                }
            
            # Normalize the query
            normalized_query = (query or "").strip().lower()
            self._samples_filters["query"] = normalized_query
            
            # Apply filters immediately
            self._apply_samples_table_filters()
            
            # Debug logging
            logging.debug(f"Search query updated: '{normalized_query}'")
        except Exception as e:
            logging.error(f"Error setting samples query: {e}")
            pass

    def _set_samples_origin_filter(self, origin_value: str) -> None:
        try:
            # Ensure we have a filters dict
            if not hasattr(self, "_samples_filters"):
                self._samples_filters = {
                    "query": "",
                    "origin": "All",
                    "job_type": "All",
                }
            
            # Normalize the origin value
            normalized_origin = origin_value or "All"
            self._samples_filters["origin"] = normalized_origin
            
            # Apply filters immediately
            self._apply_samples_table_filters()
            
            # Debug logging
            logging.debug(f"Origin filter updated: '{normalized_origin}'")
        except Exception as e:
            logging.error(f"Error setting samples origin filter: {e}")
            pass

    def _set_samples_job_type_filter(self, job_type_value: str) -> None:
        try:
            if not hasattr(self, "_samples_filters"):
                self._samples_filters = {
                    "query": "",
                    "origin": "All",
                    "job_type": "All",
                }

            normalized_job_type = job_type_value or "All"
            self._samples_filters["job_type"] = normalized_job_type

            self._apply_samples_table_filters()
            logging.debug(f"Job type filter updated: '{normalized_job_type}'")
        except Exception as e:
            logging.error(f"Error setting samples job type filter: {e}")
            pass

    def _dedupe_samples_rows_by_id(
        self, rows: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """Keep one row per sample_id (newest _last_seen_raw wins) to avoid duplicate keys."""
        by_id: Dict[str, Dict[str, Any]] = {}
        for r in rows or []:
            sid = str(r.get("sample_id", "") or "").strip()
            if not sid or sid == "unknown":
                continue
            prev = by_id.get(sid)
            if prev is None:
                by_id[sid] = r
            else:
                try:
                    if float(r.get("_last_seen_raw", 0)) >= float(
                        prev.get("_last_seen_raw", 0)
                    ):
                        by_id[sid] = r
                except Exception:
                    by_id[sid] = r
        return list(by_id.values())

    def _apply_samples_table_filters(self) -> None:
        try:
            base_rows = self._dedupe_samples_rows_by_id(
                getattr(self, "_last_samples_rows", []) or []
            )
            self._last_samples_rows = base_rows
            rows = self._normalize_rows_for_display(base_rows)
            
            logging.debug(f"Applying filters to {len(rows)} base rows")

            # Origin filter, compute dynamic 'Complete' for display if needed
            now_ts = time.time()
            for r in rows:
                    try:
                        if r.get("origin") == "Live":
                            last_raw = float(r.get("_last_seen_raw", 0))
                            active_jobs_count = r.get("active_jobs", 0)
                            pending_jobs_count = r.get("pending_jobs", 0)
                            # Only mark as Complete if timeout passed AND no active jobs
                            if last_raw and (now_ts - last_raw) >= self.completion_timeout_seconds:
                                if active_jobs_count == 0 and pending_jobs_count == 0:
                                    r["origin"] = "Complete"
                                # If there are active jobs, keep as Live even if timeout passed
                    except Exception:
                        pass

            origin = (self._samples_filters or {}).get("origin", "All")
            if origin and origin != "All":
                rows = [r for r in rows if (r.get("origin") == origin)]
                logging.debug(f"After origin filter ({origin}): {len(rows)} rows")

            # Keep job-type filter options up to date from the unfiltered row set.
            try:
                all_job_types: Set[str] = set()
                for r in self._normalize_rows_for_display(base_rows):
                    jt = str(r.get("job_types", "") or "").strip()
                    if not jt:
                        continue
                    for token in jt.split(","):
                        token = token.strip()
                        if token:
                            all_job_types.add(token)
                options = ["All"] + sorted(all_job_types)
                if hasattr(self, "job_type_filter"):
                    self.job_type_filter.set_options(options)
                    current_jt = (self._samples_filters or {}).get("job_type", "All")
                    if current_jt not in options:
                        self._samples_filters["job_type"] = "All"
                        self.job_type_filter.value = "All"
            except Exception:
                pass

            # Job type filter (exact token match in comma-separated job_types field)
            selected_job_type = (self._samples_filters or {}).get("job_type", "All")
            if selected_job_type and selected_job_type != "All":
                def _row_has_job_type(r: Dict[str, Any]) -> bool:
                    jt = str(r.get("job_types", "") or "").strip()
                    if not jt:
                        return False
                    parts = [p.strip() for p in jt.split(",") if p.strip()]
                    return selected_job_type in parts

                rows = [r for r in rows if _row_has_job_type(r)]
                logging.debug(
                    f"After job type filter ({selected_job_type}): {len(rows)} rows"
                )

            # Global query filter
            q = (self._samples_filters or {}).get("query", "")
            if q:
                ql = q.lower()
                logging.debug(f"Applying search filter: '{ql}'")

                def match_any(r: Dict[str, Any]) -> bool:
                    return any(
                        (str(r.get(k, "")).lower().find(ql) >= 0)
                        for k in ["sample_id", "origin", "job_types", "last_seen"]
                    ) or any(
                        (str(r.get(k, 0)).lower().find(ql) >= 0)
                        for k in [
                            "active_jobs",
                            "total_jobs",
                            "completed_jobs",
                            "failed_jobs",
                        ]
                    )

                rows = [r for r in rows if match_any(r)]
                logging.debug(f"After search filter: {len(rows)} rows")

            # annotate export selection state per row for rightmost checkbox column
            try:
                selected = getattr(self, "_selected_sample_ids", set()) or set()
            except Exception:
                selected = set()
            for r in rows:
                try:
                    r["export"] = r.get("sample_id") in selected
                except Exception:
                    r["export"] = False
                try:
                    sid = str(r.get("sample_id") or "")
                    s = self._sample_pipeline_status.get(sid, {})
                    r["pipeline_status"] = s.get("phase", "")
                    r["pipeline_progress"] = float(s.get("progress", 0.0) or 0.0)
                    r["pipeline_detail"] = s.get("detail", "")
                except Exception:
                    r["pipeline_status"] = ""
                    r["pipeline_progress"] = 0.0
                    r["pipeline_detail"] = ""

            # Sort rows by last activity in reverse chronological order (newest first)
            try:
                rows.sort(key=lambda r: float(r.get("_last_seen_raw", 0)), reverse=True)
            except Exception:
                # Fallback to sorting by last_seen string if _last_seen_raw is not available
                try:
                    rows.sort(key=lambda r: r.get("last_seen", ""), reverse=True)
                except Exception:
                    pass

            # Update the table
            if hasattr(self, "samples_table"):
                self.samples_table.rows = rows
                self.samples_table.update()
                logging.debug(f"Updated samples table with {len(rows)} filtered rows")
            else:
                logging.warning("samples_table not found, cannot update")
                
        except Exception as e:
            logging.error(f"Error applying samples table filters: {e}")
            pass

    def _refresh_sample_plots(self, sample_id: str):
        """Refresh all plots for a specific sample by clearing all cached state."""
        try:
            if not self.monitored_directory:
                return

            sample_dir = Path(self.monitored_directory) / sample_id
            if not sample_dir.exists():
                return

            # Clear ALL cached state to force complete regeneration
            key = str(sample_dir)

            # Completely clear coverage state for this sample
            if hasattr(self, "_coverage_state") and key in self._coverage_state:
                old_state = self._coverage_state[key]
                logging.info(
                    f"Clearing coverage state for {sample_id}: {list(old_state.keys())}"
                )
                self._coverage_state[key] = {}
                logging.info(f"Cleared all coverage state for sample {sample_id}")

            # Clear other component states if they exist
            for attr in ["_mgmt_state", "_cnv_state", "_fusion_state"]:
                if hasattr(self, attr):
                    state_dict = getattr(self, attr)
                    if key in state_dict:
                        old_state = state_dict[key]
                        logging.info(
                            f"Clearing {attr} state for {sample_id}: {list(old_state.keys())}"
                        )
                        # Clear all state for these components
                        state_dict[key] = {}
                        logging.info(f"Cleared all {attr} state for sample {sample_id}")

        except Exception as e:
            logging.debug(f"Failed to refresh sample plots for {sample_id}: {e}")

    '''
    def _refresh_all_sample_plots(self):
        """Refresh all sample plots by resetting modification time caches."""
        try:
            if not self.monitored_directory:
                return

            # Get all sample directories
            monitored_path = Path(self.monitored_directory)
            if not monitored_path.exists():
                return

            sample_dirs = [d for d in monitored_path.iterdir() if d.is_dir()]

            for sample_dir in sample_dirs:
                sample_id = sample_dir.name
                self._refresh_sample_plots(sample_id)

            ui.notify(
                f"Refreshed plots for {len(sample_dirs)} samples", type="positive"
            )
            logging.info(f"Refreshed plots for {len(sample_dirs)} samples")

        except Exception as e:
            ui.notify(f"Failed to refresh all plots: {e}", type="negative")
            logging.exception("Failed to refresh all sample plots")
    '''

    def _open_view_identifiers_modal(
        self, sample_dir: Optional[Path], sample_id: str
    ) -> None:
        """Open a modal to enter sample date of birth and view decrypted sample identifiers."""
        from cryptography.fernet import InvalidToken

        encrypted = _load_manifest_encrypted_fields(sample_dir)
        with ui.dialog().props("persistent") as dialog, ui.card().classes(
            "robin-dialog-surface w-full max-w-md p-4 md:p-5"
        ):
            ui.label("View sample identifiers").classes(
                "classification-insight-heading text-headline-small mb-2"
            )
            ui.label(
                "Enter the sample date of birth to reveal additional sample identifiers if available."
            ).classes("classification-insight-foot mb-4")
            if not encrypted:
                ui.label(
                    "No identifier manifest found for this sample, or it has no encrypted fields."
                ).classes("text-sm text-amber-700 mb-4")
                ui.button("Close", on_click=dialog.close).props(
                    "color=primary no-caps outline"
                ).classes("mt-2")
                dialog.open()
                return
            dob_input = ui.date_input(
                "Date of birth",
                value=None,
            ).classes("w-full mb-4")
            result_container = ui.column().classes("w-full mt-4 gap-2")

            def on_show() -> None:
                _v = dob_input.value
                dob_val = (_v.strftime("%Y-%m-%d") if hasattr(_v, "strftime") else (str(_v).strip() if _v else ""))
                if not dob_val:
                    ui.notify("Please enter date of birth.", type="warning")
                    return
                result_container.clear()
                with result_container:
                    try:
                        fn = _decrypt_identifier_manifest_field(
                            encrypted["first_name"], dob_val
                        )
                        ln = _decrypt_identifier_manifest_field(
                            encrypted["last_name"], dob_val
                        )
                        dob_plain = _decrypt_identifier_manifest_field(
                            encrypted["dob"], dob_val
                        )
                        ui.label("First name:").classes("font-semibold text-sm")
                        ui.label(fn or "—").classes("mb-2")
                        ui.label("Last name:").classes("font-semibold text-sm")
                        ui.label(ln or "—").classes("mb-2")
                        ui.label("Date of birth:").classes("font-semibold text-sm")
                        ui.label(dob_plain or "—").classes("mb-2")
                        if "nhs_number" in encrypted:
                            nhs_plain = _decrypt_identifier_manifest_field(
                                encrypted["nhs_number"], dob_val
                            )
                            ui.label("Hospital Number:").classes("font-semibold text-sm")
                            ui.label(nhs_plain or "—").classes("mb-2")
                    except InvalidToken:
                        ui.label("Incorrect date of birth. Please try again.").classes(
                            "text-red-600 dark:text-red-400"
                        )

            with ui.row().classes("w-full gap-2 mt-2 flex-wrap"):
                ui.button("Show identifiers", on_click=on_show).props(
                    "color=primary no-caps"
                )
                ui.button("Close", on_click=dialog.close).props("flat no-caps outline")
        dialog.open()

    def _create_sample_detail_page(self, sample_id: str):
        """Create the individual sample detail page."""
        sample_dir = (
            Path(self.monitored_directory) / sample_id
            if self.monitored_directory
            else None
        )
        test_id = _get_test_id_from_manifest(sample_dir) if sample_dir else ""

        # Check if this is a page refresh vs navigation
        # Use browser storage to detect if this is a refresh
        is_page_refresh = False
        try:
            # Check if we have a flag indicating this page was recently loaded
            if hasattr(ui, 'storage') and hasattr(ui.storage, 'browser'):
                last_load_time = ui.storage.browser.get(f'sample_{sample_id}_last_load', 0)
                current_time = time.time()
                # If last load was very recent (< 2 seconds), it's likely a page refresh
                is_page_refresh = (current_time - last_load_time) < 2.0
                
                # Update the last load time
                ui.storage.browser[f'sample_{sample_id}_last_load'] = current_time
        except Exception:
            # If storage is not available, assume it's navigation (show loading)
            is_page_refresh = False
        
        # Show loading spinner unless it's a page refresh
        show_loading = not is_page_refresh
        
        async def confirm_report_generation():
            """Show a confirmation dialog before generating the report."""
            report_types = {
                "summary": "Summary Only",
                "detailed": "Detailed",
            }
            state: Dict[str, Any] = {
                "type": "detailed",
                "export_csv": False,
                "include_sample_ids": False,
                "sample_dob": "",
            }

            with ui.dialog().props("persistent") as dialog:
                with ui.card().classes(
                    "robin-dialog-surface w-96 max-w-[95vw] p-4"
                ):
                    title_label = ui.label("Generate report").classes(
                        "classification-insight-heading text-headline-small mb-3"
                    )

                    with ui.column():
                        with ui.column().classes("mb-4"):
                            ui.label("Report type").classes(
                                "target-coverage-panel__meta-label mb-2"
                            )
                            type_toggle = ui.toggle(
                                report_types,
                                value="detailed",
                                on_change=lambda e: state.update({"type": e.value}),
                            )

                        with ui.column().classes("mb-4"):
                            ui.label("Sample identifiers").classes(
                                "target-coverage-panel__meta-label mb-2"
                            )
                            ui.label(
                                "Enter the sample date of birth to reveal additional sample identifiers if available."
                            ).classes("classification-insight-foot mb-2")

                            def on_include_sample_ids_change(e):
                                state["include_sample_ids"] = bool(e.value)
                                sample_dob_input.set_visibility(bool(e.value))

                            ui.checkbox(
                                "Include sample identifiers in report",
                                value=False,
                                on_change=on_include_sample_ids_change,
                            )
                            sample_dob_input = ui.date_input(
                                "Sample date of birth (required to reveal identifiers)",
                                value=None,
                            ).classes("w-full")
                            sample_dob_input.set_visibility(False)
                            def _on_sample_dob_change(_):
                                v = getattr(sample_dob_input, "value", None)
                                if hasattr(v, "strftime"):
                                    state["sample_dob"] = v.strftime("%Y-%m-%d")
                                else:
                                    state["sample_dob"] = str(v).strip() if v else ""
                            sample_dob_input.on("update:model-value", _on_sample_dob_change)

                        with ui.column().classes("mb-4"):
                            ui.label("Disclaimer").classes(
                                "target-coverage-panel__meta-label mb-2"
                            )
                            formatted_text = EXTENDED_DISCLAIMER_TEXT.replace(
                                "\n\n", "<br><br>"
                            ).replace("\n", " ")
                            ui.label(formatted_text).classes(
                                "text-sm text-gray-600 mb-4"
                            )

                        with ui.column().classes("mb-4"):
                            ui.label("Include data").classes(
                                "target-coverage-panel__meta-label mb-2"
                            )
                            csv_checkbox = ui.checkbox(
                                "CSV data (ZIP)",
                                value=False,
                                on_change=lambda e: state.update(
                                    {"export_csv": bool(e.value)}
                                ),
                            )

                        ui.label(
                            "Are you sure you want to generate a report?"
                        ).classes("classification-insight-foot mb-4")

                        def _capture_dob_and_confirm():
                            """Capture date picker value into state before closing dialog."""
                            if state.get("include_sample_ids"):
                                v = getattr(sample_dob_input, "value", None)
                                if v is not None:
                                    if hasattr(v, "strftime"):
                                        state["sample_dob"] = v.strftime("%Y-%m-%d")
                                    else:
                                        state["sample_dob"] = str(v).strip()
                                else:
                                    state["sample_dob"] = ""
                            dialog.submit("Yes")

                        with ui.row().classes("justify-end gap-2 flex-wrap"):
                            ui.button(
                                "No", on_click=lambda: dialog.submit("No")
                            ).props("flat no-caps outline")
                            ui.button(
                                "Yes",
                                on_click=_capture_dob_and_confirm,
                            ).props("color=primary no-caps")
                    
                    # Run name below buttons
                    ui.label(f"Exporting data for run: {sample_id}").classes(
                        "text-sm font-medium text-gray-700 mt-4"
                    )

            dialog_result = await dialog
            if dialog_result != "Yes":
                return

            # If user requested sample identifiers, decrypt with DOB
            state["sample_identifiers"] = None
            if state.get("include_sample_ids"):
                dob_val = (state.get("sample_dob") or "").strip()
                if not dob_val:
                    ui.notify(
                        "Date of birth required to include sample identifiers. Report will be generated without them.",
                        type="warning",
                    )
                else:
                    encrypted = _load_manifest_encrypted_fields(sample_dir)
                    if encrypted:
                        try:
                            from cryptography.fernet import InvalidToken
                            decrypted = {
                                "first_name": _decrypt_identifier_manifest_field(
                                    encrypted["first_name"], dob_val
                                ),
                                "last_name": _decrypt_identifier_manifest_field(
                                    encrypted["last_name"], dob_val
                                ),
                                "dob": _decrypt_identifier_manifest_field(
                                    encrypted["dob"], dob_val
                                ),
                            }
                            if "nhs_number" in encrypted:
                                decrypted["nhs_number"] = _decrypt_identifier_manifest_field(
                                    encrypted["nhs_number"], dob_val
                                )
                            decrypted["sample_id"] = sample_id
                            decrypted["test_id"] = _get_test_id_from_manifest(sample_dir)
                            state["sample_identifiers"] = decrypted
                        except InvalidToken:
                            ui.notify(
                                "Incorrect date of birth. Report will be generated without sample identifiers.",
                                type="warning",
                            )
                    else:
                        ui.notify(
                            "No identifier manifest found for this sample. Report will be generated without sample identifiers.",
                            type="info",
                        )

            # Now show the progress dialog
            with ui.dialog().props("persistent") as progress_dialog:
                with ui.card().classes(
                    "robin-dialog-surface w-96 max-w-[95vw] p-4"
                ):
                    ui.label("Generating report").classes(
                        "classification-insight-heading text-headline-small mb-2"
                    )
                    
                    # Divider
                    ui.separator().classes("mb-4")
                    
                    # Run name
                    ui.label(f"Exporting data for run: {sample_id}").classes(
                        "text-sm font-medium text-gray-700 mb-4"
                    )

                    # Report type and output displays
                    report_type_display = ui.label(
                        f"Report Type: {state.get('type', 'detailed').title()}"
                    ).classes("text-sm mb-2")
                    
                    output_display = ui.label(
                        f"Output: PDF{' + CSV (ZIP)' if state.get('export_csv', False) else ''}"
                    ).classes("text-sm mb-4")
                    
                    # Progress indicator
                    progress_bar = ui.linear_progress(0.0).classes("mb-2")
                    
                    # Progress text
                    progress_text = ui.label("Preparing...").classes(
                        "text-xs text-gray-600 text-center mb-2"
                    )
                    
                    # Messages container - use label with newlines for multiple messages
                    messages_label = ui.label("").classes("text-xs text-gray-500")

                    # Track messages
                    progress_updates = queue.Queue()
                    messages_list = []
                    
                    # Timer to process progress updates on UI thread
                    def process_progress_updates():
                        """Process queued progress updates."""
                        try:
                            while not progress_updates.empty():
                                update = progress_updates.get_nowait()
                                stage = update.get("stage", "unknown")
                                message = update.get("message", "")
                                progress = update.get("progress", 0.0)
                                
                                if progress is not None:
                                    progress_bar.value = progress
                                    progress_text.text = f"{int(progress * 100)}% - {message}"
                                else:
                                    progress_text.text = message
                                
                                # Add message to messages list
                                messages_list.append(message)
                                # Keep only last 10 messages
                                if len(messages_list) > 10:
                                    messages_list.pop(0)
                                
                                # Update messages label
                                messages_label.text = "\n".join(messages_list[-5:])
                                
                        except queue.Empty:
                            pass
                        except Exception as e:
                            logging.debug(f"Error processing progress updates: {e}")
                    
                    # Set up timer to process updates
                    update_timer = ui.timer(0.1, process_progress_updates)
                    
                    def progress_callback(progress_data: Dict[str, Any]):
                        """Custom progress callback to update dialog (called from background thread)."""
                        try:
                            # Queue the update instead of directly updating UI
                            progress_updates.put(progress_data)
                        except Exception as e:
                            logging.debug(f"Error in progress callback: {e}")
                    
                    # Track if still generating
                    is_generating = {"active": True}

                    # Storage for the files to download (will be populated by background task)
                    files_to_download = []
                    download_complete = {"done": False}
                    
                    # Timer to handle downloads once background task is done
                    def handle_downloads():
                        """Handle downloads in UI context once generation is complete."""
                        if download_complete["done"] and files_to_download:
                            # Log that report generation is complete and downloads are available
                            file_count = len([f for f in files_to_download if f is not None])
                            logging.info(f"Report generation complete for {sample_id}. {file_count} file(s) ready for download.")
                            
                            for file_path in files_to_download:
                                if file_path is not None:
                                    logging.debug(f"Initiating download: {file_path}")
                                    ui.download(file_path)
                            # Close dialog after 3 seconds
                            ui.timer(3.0, lambda: progress_dialog.submit(None), once=True)
                            download_timer.deactivate()
                    
                    download_timer = ui.timer(0.1, handle_downloads)
                    
                    # Start report generation
                    async def complete_generation():
                        """Complete the report generation and close dialog."""
                        try:
                            await generate_and_download_report(state, progress_callback, progress_dialog, is_generating, files_to_download)
                        finally:
                            # Clean up timer
                            update_timer.deactivate()
                            is_generating["active"] = False
                            download_complete["done"] = True
                    
                    asyncio.create_task(complete_generation())

            await progress_dialog

        async def generate_and_download_report(state: Dict[str, Any], progress_callback, progress_dialog, is_generating, files_to_download):
            """Generate report and update progress in dialog."""
            try:
                from nicegui import run as ng_run  # type: ignore
                
                if not sample_dir or not sample_dir.exists():
                    # Queue error notification
                    files_to_download.append(None)  # Signal error
                    ui.timer(0.1, lambda: ui.notify(
                        "Output directory not available for this sample",
                        type="warning",
                    ), once=True)
                    return
                
                filename = f"{sample_id}_run_report.pdf"
                pdf_path = os.path.join(str(sample_dir), filename)
                os.makedirs(str(sample_dir), exist_ok=True)
                export_csv_dir = None
                if bool(state.get("export_csv", False)):
                    export_csv_dir = os.path.join(
                        str(sample_dir), "report_csv"
                    )
                
                # Use only our custom callback, not the notification system
                def combined_callback(progress_data: Dict[str, Any]):
                    """Custom callback for dialog updates only (no notifications)."""
                    progress_callback(progress_data)
                
                pdf_file = await ng_run.io_bound(
                    create_pdf,
                    pdf_path,
                    str(sample_dir),
                    self.center or "Unknown",
                    report_type=state.get("type", "detailed"),
                    export_csv_dir=export_csv_dir,
                    export_xlsx=False,
                    export_zip=bool(state.get("export_csv", False)),
                    progress_callback=combined_callback,
                    workflow_steps=self.workflow_steps if hasattr(self, 'workflow_steps') else None,
                    sample_identifiers=state.get("sample_identifiers"),
                )
                
                # Mark report as completed
                from robin.gui.report_progress import progress_manager
                progress_manager.complete_report(sample_id, filename)
                
                # Queue files for download in UI context
                files_to_download.append(pdf_file)
                
                # Also offer CSV ZIP if requested
                if bool(state.get("export_csv", False)) and export_csv_dir:
                    zip_path = os.path.join(
                        export_csv_dir, f"{sample_id}_report_data.zip"
                    )
                    if os.path.exists(zip_path):
                        files_to_download.append(zip_path)
                
            except Exception as e:
                # Mark report as failed
                from robin.gui.report_progress import progress_manager
                progress_manager.error_report(sample_id, str(e))
                
                # Queue error notification in UI context
                ui.timer(0.1, lambda: ui.notify(
                    f"Error generating report: {str(e)}",
                    type="negative",
                ), once=True)
                files_to_download.append(None)  # Signal error

        async def download_report(state: Dict[str, Any]):
            """Generate and download the report for this sample."""
            try:
                # Import here to avoid global dependency if GUI isn't used
                from nicegui import run as ng_run  # type: ignore

                
                if not sample_dir or not sample_dir.exists():
                    ui.notify(
                        "Output directory not available for this sample",
                        type="warning",
                    )
                    return
                
                # Show initial notification with more explicit styling
                ui.notify(
                    f"[{sample_id}] Starting report generation...",
                    type="info",
                    timeout=0,  # Persistent notification
                    position="top-right"
                )
                
                filename = f"{sample_id}_run_report.pdf"
                pdf_path = os.path.join(str(sample_dir), filename)
                os.makedirs(str(sample_dir), exist_ok=True)
                export_csv_dir = None
                if bool(state.get("export_csv", False)):
                    export_csv_dir = os.path.join(
                        str(sample_dir), "report_csv"
                    )
                
                # Create progress callback
                from robin.gui.report_progress import create_progress_callback
                progress_callback = create_progress_callback(sample_id)
                
                pdf_file = await ng_run.io_bound(
                    create_pdf,
                    pdf_path,
                    str(sample_dir),
                    self.center or "Unknown",
                    report_type=state.get("type", "detailed"),
                    export_csv_dir=export_csv_dir,
                    export_xlsx=False,
                    export_zip=bool(state.get("export_csv", False)),
                    progress_callback=progress_callback,
                    workflow_steps=self.workflow_steps if hasattr(self, 'workflow_steps') else None,
                )
                
                # Mark report as completed
                from robin.gui.report_progress import progress_manager
                progress_manager.complete_report(sample_id, filename)
                
                ui.download(pdf_file)
                # Also offer CSV ZIP if requested
                if bool(state.get("export_csv", False)) and export_csv_dir:
                    zip_path = os.path.join(
                        export_csv_dir, f"{sample_id}_report_data.zip"
                    )
                    if os.path.exists(zip_path):
                        ui.download(zip_path)
                
            except Exception as e:
                # Mark report as failed
                from robin.gui.report_progress import progress_manager
                progress_manager.error_report(sample_id, str(e))
                
        title_suffix = f" | {test_id}" if test_id else ""
        with theme.frame(
            f"R.O.B.I.N - Sample {sample_id}{title_suffix}",
            smalltitle=f"{sample_id}{title_suffix}".strip(),
            batphone=False,
            center=self.center,
            setup_notifications=self._setup_notification_system,
        ):
            # Guard: unknown sample -> show message and back button; also redirect
            if self._known_sample_ids and sample_id not in self._known_sample_ids:
                with ui.column().classes(
                    "w-full max-w-lg mx-auto items-center justify-center p-5 gap-2"
                ):
                    ui.label(f"Unknown sample").classes(
                        "text-headline-medium text-slate-900 dark:text-slate-50"
                    )
                    ui.label(sample_id).classes(
                        "text-display-small text-red-600 dark:text-rose-400 font-mono"
                    )
                    ui.label(
                        "This library ID has not been seen yet in the current session."
                    ).classes("text-body-medium text-slate-600 dark:text-slate-400 text-center")
                    ui.label(
                        "Redirecting to sample list…"
                    ).classes("text-body-small text-slate-500 dark:text-slate-500 mt-1")
                    ui.button(
                        "Back to samples", on_click=lambda: ui.navigate.to("/live_data")
                    ).props("color=primary").classes("rounded-lg mt-2")
                
                # Soft redirect after short delay
                def redirect_to_samples():
                    """Redirect to the sample list page."""
                    try:
                        ui.navigate.to("/live_data")
                    except Exception as e:
                        logging.warning(f"Failed to redirect to /live_data: {e}")
                        # Fallback: try using JavaScript redirect
                        try:
                            ui.run_javascript('window.location.href = "/live_data"')
                        except Exception as e2:
                            logging.error(f"Failed JavaScript redirect: {e2}")
                
                # Create and activate the timer
                redirect_timer = ui.timer(2.0, redirect_to_samples, once=True)
                redirect_timer.activate()
                return

            # Debug visibility: entering sample page and directory availability
            try:
                ui.notify(f"Opening sample {sample_id}", type="info")
            except Exception:
                pass
            
            # Check directory existence asynchronously to avoid blocking
            async def check_directory_and_notify():
                try:
                    if not sample_dir or not sample_dir.exists():
                        ui.notify(
                            f"Sample output directory not found for {sample_id}",
                            type="warning",
                        )
                except Exception:
                    pass
            
            # Start directory check in background
            ui.timer(0.1, check_directory_and_notify, once=True)

            with ui.card().classes(
                "w-full border border-slate-200 dark:border-slate-800 rounded-xl overflow-hidden shadow-sm"
            ):
                with ui.column().classes(
                    "w-full min-w-0 gap-0 bg-gradient-to-b from-slate-50 to-white "
                    "dark:from-slate-900/90 dark:to-[var(--md-surface)] "
                    "border-b border-slate-200/90 dark:border-slate-800"
                ):
                    # Mobile: flex-col stacks title + actions; md+ row (Quasar ui.column ignores flex-row)
                    with ui.element("div").classes(
                        "w-full min-w-0 flex flex-col gap-2 p-2 md:p-3 "
                        "md:flex-row md:justify-between md:items-start"
                    ):
                        with ui.column().classes("gap-3 min-w-0 w-full md:flex-1"):
                            # Anchor: Manrope ExtraBold — break-words, not break-all (mid-word splits)
                            ui.label(sample_id).classes(
                                "text-display-small font-extrabold text-slate-900 "
                                "dark:text-slate-50 break-words w-full min-w-0"
                            )
                            # Data labels: Inter semibold; values: JetBrains Mono (§2.2)
                            with ui.column().classes("gap-1.5"):
                                with ui.row().classes("items-baseline gap-2 flex-wrap"):
                                    ui.label("Library ID").classes(
                                        "text-label-medium text-slate-500 dark:text-slate-400"
                                    )
                                    ui.label(sample_id).classes(
                                        "text-body-small text-mono text-slate-800 "
                                        "dark:text-slate-200 break-words"
                                    )
                                if test_id:
                                    with ui.row().classes("items-baseline gap-2 flex-wrap"):
                                        ui.label("Test ID").classes(
                                            "text-label-medium text-slate-500 dark:text-slate-400"
                                        )
                                        ui.label(test_id).classes(
                                            "text-body-small text-mono text-slate-800 "
                                            "dark:text-slate-200 break-words"
                                        )
                            ui.label(
                                "Classification, coverage, and analysis modules below."
                            ).classes(
                                "text-body-medium text-slate-600 dark:text-slate-400"
                            )
                            ui.button(
                                "View sample identifiers",
                                on_click=lambda: self._open_view_identifiers_modal(
                                    sample_dir, sample_id
                                ),
                            ).classes(
                                "mt-1 w-full sm:w-auto self-stretch sm:self-start"
                            ).props(
                                "flat dense color=primary no-caps"
                            )
                        target_bam_exists = False
                        if sample_dir and sample_dir.exists():
                            target_bam = sample_dir / "target.bam"
                            target_bam_exists = (
                                target_bam.exists() and target_bam.is_file()
                            )
                        with ui.column().classes(
                            "w-full shrink-0 flex flex-col gap-2 "
                            "md:w-auto md:min-w-[12rem] md:items-end"
                        ):
                            if target_bam_exists:
                                ui.button(
                                    "More details",
                                    on_click=lambda: ui.navigate.to(
                                        f"/live_data/{sample_id}/details"
                                    ),
                                ).classes(
                                    "rounded-lg border border-slate-300 dark:border-slate-600 "
                                    "text-title-medium w-full md:w-auto md:min-w-[10rem]"
                                ).props("flat no-caps")
                            ui.button(
                                "Generate report",
                                on_click=confirm_report_generation,
                            ).props("color=primary no-caps").classes(
                                "rounded-lg px-4 py-2 text-title-medium "
                                "w-full md:w-auto md:min-w-[10rem]"
                            )
                # Main content area with conditional loading state
                with ui.column().classes("w-full p-2 md:p-3 gap-2"):
                    if show_loading:
                        # Loading container that will be hidden when data is ready
                        loading_container = ui.column().classes(
                            "w-full items-center justify-center p-5 "
                            "rounded-xl border border-slate-200/90 dark:border-slate-800 "
                            "bg-gradient-to-b from-slate-50 to-white "
                            "dark:from-slate-900/60 dark:to-zinc-950/80"
                        )
                        with loading_container:
                            ui.spinner("bars", size="4em", color="primary").classes("mb-4")
                            ui.label("Loading sample data…").classes(
                                "text-title-medium text-slate-800 dark:text-slate-100"
                            )
                            ui.label("This may take a moment for large datasets").classes(
                                "text-body-small text-slate-500 dark:text-slate-400"
                            )

                        # Content container that will be shown when data is ready
                        content_container = (
                            ui.column()
                            .classes("w-full gap-2")
                            .style("display: none")
                        )
                    else:
                        # For page refreshes, show content immediately
                        content_container = ui.column().classes("w-full gap-2")
                        loading_container = None

                    with content_container:
                        # Summary section (new component) - create UI immediately, load data async
                        try:
                            try:
                                from .gui.components.summary import add_summary_section  # type: ignore
                            except ImportError:
                                # Try absolute import if relative fails
                                from robin.gui.components.summary import add_summary_section

                            # Create the UI components immediately on the main thread
                            add_summary_section(sample_dir, sample_id, self)
                        except Exception as e:
                            logging.exception(f"[GUI] Summary section failed: {e}")
                            try:
                                ui.notify(f"Summary section failed: {e}", type="warning")
                            except Exception:
                                pass

                        # Import section visibility helper
                        try:
                            from robin.gui.config import is_section_enabled, get_enabled_classification_steps
                        except ImportError:
                            is_section_enabled = lambda name, steps: True
                            get_enabled_classification_steps = lambda steps: {"sturgeon", "nanodx", "random_forest", "pannanodx"}
                        
                        workflow_steps = self.workflow_steps if hasattr(self, 'workflow_steps') else None

                        # Defer heavy analysis sections until after initial render
                        analysis_loading_card = ui.card().classes(
                            "w-full p-2 border border-slate-200/90 dark:border-slate-800 "
                            "bg-gradient-to-b from-white to-slate-50/80 "
                            "dark:from-slate-900/40 dark:to-zinc-950/60"
                        )
                        with analysis_loading_card:
                            with ui.row().classes("items-center gap-3"):
                                ui.spinner(size="md", color="primary")
                                ui.label("Loading analysis sections…").classes(
                                    "text-body-medium text-slate-700 dark:text-slate-200"
                                )
                            ui.label(
                                "This can take a moment for large datasets."
                            ).classes("text-body-small text-slate-500 dark:text-slate-400 mt-2")

                        analysis_container = ui.column().classes("w-full gap-2")

                        def _build_analysis_sections():
                            try:
                                try:
                                    analysis_loading_card.set_visibility(False)
                                except Exception:
                                    try:
                                        analysis_loading_card.style("display: none")
                                    except Exception:
                                        pass
                                analysis_container.clear()
                                with analysis_container:
                                    # MNP-Flex section
                                    try:
                                        try:
                                            from .gui.components.mnpflex import add_mnpflex_section  # type: ignore
                                        except ImportError:
                                            from robin.gui.components.mnpflex import add_mnpflex_section

                                        add_mnpflex_section(self, sample_dir, sample_id)
                                    except Exception as e:
                                        logging.exception(f"[GUI] MNP-Flex section failed: {e}")
                                        try:
                                            ui.notify(
                                                f"MNP-Flex section failed: {e}", type="warning"
                                            )
                                        except Exception:
                                            pass

                                    # Classification section
                                    enabled_classification_steps = get_enabled_classification_steps(workflow_steps)
                                    if not workflow_steps or enabled_classification_steps:
                                        try:
                                            try:
                                                from .gui.components.classification import add_classification_section  # type: ignore
                                            except ImportError:
                                                from robin.gui.components.classification import (
                                                    add_classification_section,
                                                )

                                            add_classification_section(sample_dir, self)
                                        except Exception as e:
                                            logging.exception(f"[GUI] Classification section failed: {e}")
                                            try:
                                                ui.notify(
                                                    f"Classification section failed: {e}", type="warning"
                                                )
                                            except Exception:
                                                pass

                                    # Coverage section (target)
                                    if not workflow_steps or is_section_enabled("target", workflow_steps):
                                        try:
                                            try:
                                                from .gui.components.coverage import add_coverage_section  # type: ignore
                                            except ImportError:
                                                from robin.gui.components.coverage import (
                                                    add_coverage_section,
                                                    add_igv_viewer,
                                                )

                                            add_coverage_section(self, sample_dir)
                                        except Exception as e:
                                            try:
                                                ui.notify(f"Coverage section failed: {e}", type="warning")
                                            except Exception:
                                                pass
                                    
                                    # MGMT section
                                    if not workflow_steps or is_section_enabled("mgmt", workflow_steps):
                                        try:
                                            try:
                                                from .gui.components.mgmt import add_mgmt_section  # type: ignore
                                            except ImportError:
                                                from robin.gui.components.mgmt import add_mgmt_section

                                            add_mgmt_section(self, sample_dir)
                                        except Exception as e:
                                            logging.exception(f"[GUI] MGMT section failed: {e}")
                                            try:
                                                ui.notify(f"MGMT section failed: {e}", type="warning")
                                            except Exception:
                                                pass
                                    
                                    # CNV section
                                    if not workflow_steps or is_section_enabled("cnv", workflow_steps):
                                        try:
                                            try:
                                                from .gui.components.cnv import add_cnv_section  # type: ignore
                                            except ImportError:
                                                from robin.gui.components.cnv import add_cnv_section

                                            add_cnv_section(self, sample_dir)
                                        except Exception as e:
                                            logging.exception(f"[GUI] CNV section failed: {e}")
                                            try:
                                                ui.notify(f"CNV section failed: {e}", type="warning")
                                            except Exception:
                                                pass

                                    # Fusion section + BED coverage
                                    if not workflow_steps or is_section_enabled("fusion", workflow_steps):
                                        try:
                                            try:
                                                from .gui.components.fusion import add_fusion_section  # type: ignore
                                            except ImportError:
                                                from robin.gui.components.fusion import add_fusion_section

                                            add_fusion_section(self, sample_dir)
                                        except Exception as e:
                                            logging.exception(f"[GUI] Fusion section failed: {e}")
                                            try:
                                                ui.notify(f"Fusion section failed: {e}", type="warning")
                                            except Exception:
                                                pass
                                        
                                        try:
                                            try:
                                                from .gui.components.bed_coverage import add_bed_coverage_section  # type: ignore
                                            except ImportError:
                                                from robin.gui.components.bed_coverage import add_bed_coverage_section
                                            
                                            add_bed_coverage_section(self, sample_dir)
                                        except Exception as e:
                                            logging.exception(f"[GUI] BED Coverage section failed: {e}")
                                            try:
                                                ui.notify(f"BED Coverage section failed: {e}", type="warning")
                                            except Exception:
                                                pass
                            except Exception as e:
                                logging.exception(f"[GUI] Failed to build analysis sections: {e}")

                        # Delay heavy UI creation to allow websocket handshake to complete
                        ui.timer(0.5, _build_analysis_sections, once=True)

                        from robin.gui.theme import styled_table

                        # Files in output directory (design.md §9 — insight shell)
                        with ui.element("div").classes("w-full min-w-0").props(
                            "id=analysis-detail-output-files"
                        ):
                            with ui.element("div").classes(
                                "classification-insight-shell w-full min-w-0"
                            ):
                                ui.label("Output files").classes(
                                    "classification-insight-heading text-headline-small"
                                )
                                with ui.element("div").classes(
                                    "classification-insight-card w-full min-w-0"
                                ):
                                    with ui.column().classes(
                                        "w-full min-w-0 gap-2 p-2 md:p-3"
                                    ):
                                        with ui.row().classes(
                                            "items-center gap-2 min-w-0"
                                        ):
                                            ui.icon("folder_open").classes(
                                                "classification-insight-icon"
                                            )
                                            ui.label("Sample output directory").classes(
                                                "classification-insight-model flex-1 min-w-0"
                                            )
                                        ui.label(sample_id).classes(
                                            "classification-insight-result w-full"
                                        )
                                        ui.label(str(sample_dir)).classes(
                                            "classification-insight-meta w-full break-all"
                                        )
                                        ui.label(
                                            "Browse or download files from this run."
                                        ).classes("classification-insight-foot")

                                ui.label("File list").classes(
                                    "target-coverage-panel__meta-label mt-4 mb-1"
                                )
                                with ui.row().classes(
                                    "items-center gap-3 mb-2 w-full flex-wrap"
                                ):
                                    files_search = ui.input("Search files…").props(
                                        "borderless dense clearable"
                                    )

                                _files_container, files_table = styled_table(
                                    columns=[
                                        {
                                            "name": "name",
                                            "label": "File",
                                            "field": "name",
                                            "sortable": True,
                                        },
                                        {
                                            "name": "size",
                                            "label": "Size (bytes)",
                                            "field": "size",
                                            "sortable": True,
                                        },
                                        {
                                            "name": "mtime",
                                            "label": "Last Modified",
                                            "field": "mtime",
                                            "sortable": True,
                                        },
                                        {
                                            "name": "actions",
                                            "label": "Download",
                                            "field": "actions",
                                            "sortable": False,
                                            "align": "center",
                                        },
                                    ],
                                    rows=[],
                                    pagination=20,
                                    class_size="table-xs",
                                )
                                try:
                                    files_table.props(
                                        "multi-sort rows-per-page-options="
                                        '"[10,20,50,0]"'
                                    )
                                    files_search.bind_value(files_table, "filter")
                                except Exception:
                                    pass

                                # Add download button slot that emits an event
                                try:
                                    files_table.add_slot(
                                        "body-cell-actions",
                                        """
<q-td key="actions" :props="props">
  <q-btn color="primary" size="sm" icon="download"
         @click="() => $parent.$emit('download-file', props.row.name)"
         title="Download file" />
</q-td>
""",
                                    )
                                except Exception:
                                    pass

                                try:
                                    files_table.on(
                                        "download-file",
                                        lambda event: _download_file(event.args),
                                    )
                                except Exception:
                                    pass

                    # master.csv fields (BAM counters, etc.) are shown in Run summary (summary.py)

                    # Download method for individual files
                    def _download_file(filename: str):
                        """Download a single file from the sample directory."""
                        try:
                            if not sample_dir or not sample_dir.exists():
                                ui.notify("Sample directory not found", type="error")
                                return

                            file_path = sample_dir / filename
                            if not file_path.exists() or not file_path.is_file():
                                ui.notify(f"File {filename} not found", type="error")
                                return

                            with open(file_path, 'rb') as f:
                                content = f.read()

                            ui.download(
                                content,
                                filename=filename,
                                media_type='application/octet-stream'
                            )

                        except Exception as e:
                            ui.notify(f"Download failed: {e}", type="error")

                    # Periodic refresher for files table
                    _notify_state = {"files_error": False}

                    def _refresh_files_list_sync() -> List[Dict[str, Any]]:
                        """Synchronous file list refresh - runs in background thread"""
                        rows = []
                        if sample_dir and sample_dir.exists():
                            for f in sorted(sample_dir.iterdir()):
                                if f.is_file():
                                    try:
                                        stat = f.stat()
                                        rows.append(
                                            {
                                                "name": f.name,
                                                "size": stat.st_size,
                                                "mtime": time.strftime(
                                                    "%Y-%m-%d %H:%M:%S",
                                                    time.localtime(stat.st_mtime),
                                                ),
                                                "actions": f.name,  # Store filename for actions
                                            }
                                        )
                                    except Exception:
                                        continue
                        return rows

                    async def _refresh_sample_detail_async() -> None:
                        """Asynchronous version of sample detail refresh"""
                        try:
                            # Run file operations in background threads
                            #import concurrent.futures
                            # Run file I/O operations in background threads to avoid blocking GUI
                            import asyncio
                            files_result = await asyncio.to_thread(_refresh_files_list_sync)

                            # Update UI with results
                            files_table.rows = files_result
                            files_table.update()

                            # Reset error states on success
                            _notify_state["files_error"] = False
                            
                        except Exception as e:
                            logging.error(f"Error in async sample detail refresh: {e}")
                            # Only show error notification once per error type
                            if not _notify_state["files_error"]:
                                try:
                                    ui.notify(
                                        f"Failed to refresh sample data for {sample_id}: {e}",
                                        type="warning",
                                    )
                                except Exception:
                                    pass
                                _notify_state["files_error"] = True

                    def _refresh_sample_detail() -> None:
                        """Synchronous version - kept for backward compatibility"""
                        # Refresh files list
                        try:
                            rows = _refresh_files_list_sync()
                            files_table.rows = rows
                            files_table.update()
                        except Exception as e:
                            if not _notify_state["files_error"]:
                                try:
                                    ui.notify(
                                        f"Failed to list output files for {sample_id}: {e}",
                                        type="warning",
                                    )
                                except Exception:
                                    pass
                                _notify_state["files_error"] = True

                    # Show content and hide loading after initial data load (only when showing loading)
                    if show_loading and loading_container:
                        async def _load_initial_data_and_show():
                            """Load initial data asynchronously then show content"""
                            try:
                                # Load initial data
                                await _refresh_sample_detail_async()
                                
                                # Show content and hide loading
                                loading_container.style("display: none")
                                content_container.style("display: flex")
                                
                            except Exception as e:
                                logging.error(f"Error loading initial data: {e}")
                                # Show content anyway to avoid infinite loading
                                loading_container.style("display: none")
                                content_container.style("display: flex")

                        # Start initial data loading
                        try:
                            ui.timer(0.1, _load_initial_data_and_show, once=True)
                        except Exception:
                            # Fallback: show content immediately if timer fails
                            loading_container.style("display: none")
                            content_container.style("display: flex")
                    else:
                        # For page refreshes, load data immediately
                        try:
                            ui.timer(0.1, _refresh_sample_detail_async, once=True)
                        except Exception:
                            pass

                    # Start periodic refresh with async version
                    try:
                        refresh_timer = ui.timer(30.0, _refresh_sample_detail_async)
                        try:
                            ui.context.client.on_disconnect(
                                lambda: refresh_timer.deactivate()
                            )
                        except Exception:
                            pass
                    except Exception:
                        pass

    def _create_sample_details_page(self, sample_id: str):
        """Create the sample details page with comprehensive information."""
        # Get sample directory and test_id from manifest if present
        sample_dir = (
            Path(self.monitored_directory) / sample_id
            if self.monitored_directory
            else None
        )
        test_id = _get_test_id_from_manifest(sample_dir) if sample_dir else ""

        # Check if sample is known
        if self._known_sample_ids and sample_id not in self._known_sample_ids:
            with theme.frame(
                f"R.O.B.I.N - Sample Details",
                smalltitle="Details",
                batphone=False,
                center=self.center,
                setup_notifications=self._setup_notification_system,
            ):
                with ui.column().classes(
                    "w-full max-w-lg mx-auto items-center justify-center p-5 gap-2"
                ):
                    ui.label("Unknown sample").classes(
                        "text-headline-medium text-slate-900 dark:text-slate-50"
                    )
                    ui.label(sample_id).classes(
                        "text-display-small text-red-600 dark:text-rose-400 font-mono"
                    )
                    ui.label(
                        "This library ID has not been seen yet in the current session."
                    ).classes(
                        "text-body-medium text-slate-600 dark:text-slate-400 text-center"
                    )
                    ui.button(
                        "Back to sample",
                        on_click=lambda: ui.navigate.to(f"/live_data/{sample_id}"),
                    ).props("color=primary").classes("rounded-lg mt-2")
                    ui.button(
                        "Back to samples",
                        on_click=lambda: ui.navigate.to("/live_data"),
                    ).classes(
                        "mt-1 rounded-lg border border-slate-300 dark:border-slate-600"
                    ).props("flat")
            return
        
        # Create the page with theme frame
        title_suffix = f" | {test_id}" if test_id else ""
        with theme.frame(
            f"R.O.B.I.N - Sample Details: {sample_id}{title_suffix}",
            smalltitle=f"{sample_id} Details{title_suffix}".strip(),
            batphone=False,
            center=self.center,
            setup_notifications=self._setup_notification_system,
        ):
            with ui.element("div").classes("w-full min-w-0").props("id=sample-details-page"):
                with ui.element("div").classes("classification-insight-shell w-full min-w-0"):
                    with ui.row().classes(
                        "w-full flex flex-col gap-3 md:flex-row md:justify-between "
                        "md:items-start p-2 md:p-3"
                    ):
                        with ui.column().classes("gap-2 min-w-0 flex-1"):
                            ui.label("Sample details").classes(
                                "classification-insight-heading text-headline-small"
                            )
                            ui.label(sample_id).classes(
                                "classification-insight-result w-full break-words"
                            )
                            details_label = f"Library ID · {sample_id}"
                            if test_id:
                                details_label += f" · Test ID {test_id}"
                            ui.label(details_label).classes(
                                "classification-insight-meta w-full font-mono break-all"
                            )
                            ui.label(
                                "IGV browser, sample identifiers, SNP tables, fusion pairs, "
                                "and target genes."
                            ).classes("classification-insight-foot")
                            ui.button(
                                "View sample identifiers",
                                on_click=lambda: self._open_view_identifiers_modal(
                                    sample_dir, sample_id
                                ),
                            ).classes("mt-1 self-start").props(
                                "flat dense color=primary no-caps"
                            )
                        ui.button(
                            "Back to sample",
                            on_click=lambda: ui.navigate.to(f"/live_data/{sample_id}"),
                        ).props("color=primary no-caps").classes(
                            "rounded-lg shrink-0 self-stretch md:self-start "
                            "w-full md:w-auto md:min-w-[10rem]"
                        )

                with ui.column().classes("w-full gap-3 p-2 md:p-3"):
                    with ui.element("div").classes(
                        "classification-insight-shell w-full min-w-0"
                    ):
                        ui.label("Output location").classes(
                            "classification-insight-heading text-headline-small"
                        )
                        with ui.element("div").classes(
                            "classification-insight-card w-full min-w-0"
                        ):
                            with ui.column().classes(
                                "w-full min-w-0 gap-2 p-2 md:p-3"
                            ):
                                with ui.row().classes("items-center gap-2 min-w-0"):
                                    ui.icon("folder_open").classes(
                                        "classification-insight-icon"
                                    )
                                    ui.label("Sample output directory").classes(
                                        "classification-insight-model flex-1 min-w-0"
                                    )
                                if sample_dir and sample_dir.exists():
                                    ui.label(str(sample_dir)).classes(
                                        "classification-insight-meta w-full font-mono break-all"
                                    )
                                    ui.label("Directory found").classes(
                                        "classification-insight-level "
                                        "classification-insight-level--high w-full"
                                    )
                                else:
                                    ui.label("Directory not found").classes(
                                        "classification-insight-level "
                                        "classification-insight-level--low w-full"
                                    )
                                    if sample_dir:
                                        ui.label(
                                            f"Expected path: {sample_dir}"
                                        ).classes("classification-insight-foot")

                    if self.center:
                        with ui.element("div").classes(
                            "classification-insight-shell w-full min-w-0"
                        ):
                            ui.label("Analysis center").classes(
                                "classification-insight-heading text-headline-small"
                            )
                            with ui.element("div").classes(
                                "classification-insight-card w-full min-w-0"
                            ):
                                with ui.column().classes(
                                    "w-full min-w-0 gap-2 p-2 md:p-3"
                                ):
                                    with ui.row().classes(
                                        "items-center gap-2 min-w-0"
                                    ):
                                        ui.icon("hub").classes(
                                            "classification-insight-icon"
                                        )
                                        ui.label("Deployment").classes(
                                            "classification-insight-model flex-1 min-w-0"
                                        )
                                    ui.label(str(self.center)).classes(
                                        "classification-insight-result w-full break-words"
                                    )

                    # IGV Viewer section - moved to top, before tables
                    if sample_dir and sample_dir.exists():
                        from robin.gui.components.coverage import add_igv_viewer
                        add_igv_viewer(self, sample_dir)
                    
                    # SNP Analysis section
                    if sample_dir and sample_dir.exists():
                        from robin.gui.components.snp import add_snp_section
                        add_snp_section(self, sample_dir)
                    
                    # Fusion Pairs Table section
                    if sample_dir and sample_dir.exists():
                        from robin.gui.components.fusion import (
                            _load_processed_pickle,
                            _cluster_fusion_reads
                        )
                        import pandas as pd
                        
                        # Load fusion data
                        fusion_data_loaded = False
                        fusion_data = None
                        try:
                            target_file = sample_dir / "fusion_candidates_master_processed.pkl"
                            genome_file = sample_dir / "fusion_candidates_all_processed.pkl"
                            
                            # Try to load target panel data first, then genome-wide
                            if target_file.exists():
                                fusion_data = _load_processed_pickle(target_file)
                                if fusion_data and fusion_data.get("annotated_data") is not None:
                                    fusion_data_loaded = True
                            elif genome_file.exists():
                                fusion_data = _load_processed_pickle(genome_file)
                                if fusion_data and fusion_data.get("annotated_data") is not None:
                                    fusion_data_loaded = True
                        except Exception as e:
                            logging.warning(f"Failed to load fusion data: {e}")
                        
                        if fusion_data_loaded and fusion_data:
                            annotated_data = fusion_data.get("annotated_data", pd.DataFrame())
                            goodpairs = fusion_data.get("goodpairs", pd.Series())
                            
                            if not annotated_data.empty:
                                # Filter to good pairs if available
                                if not goodpairs.empty and goodpairs.sum() > 0:
                                    aligned_goodpairs = goodpairs.reindex(annotated_data.index, fill_value=False)
                                    filtered_data = annotated_data[aligned_goodpairs]
                                else:
                                    filtered_data = annotated_data
                                
                                # Cluster fusion reads
                                clustered_data = _cluster_fusion_reads(
                                    filtered_data, 
                                    max_distance=10000, 
                                    use_breakpoint_validation=True
                                )
                                
                                if not clustered_data.empty:
                                    with ui.element("div").classes(
                                        "classification-insight-shell w-full min-w-0"
                                    ):
                                        ui.label("Fusion pairs").classes(
                                            "classification-insight-heading text-headline-small"
                                        )
                                        ui.label(
                                            "Click a row to open the fusion region in IGV."
                                        ).classes("classification-insight-meta w-full mb-2")
                                        
                                        # Create columns for the table
                                        from robin.gui.theme import styled_table
                                        
                                        columns = [
                                            {"name": "fusion_pair", "label": "Fusion Pair", "field": "fusion_pair", "sortable": True},
                                            {"name": "chr1", "label": "Chr 1", "field": "chr1", "sortable": True},
                                            {"name": "pos1", "label": "Breakpoint 1", "field": "pos1", "sortable": True},
                                            {"name": "chr2", "label": "Chr 2", "field": "chr2", "sortable": True},
                                            {"name": "pos2", "label": "Breakpoint 2", "field": "pos2", "sortable": True},
                                            {"name": "reads", "label": "Supporting Reads", "field": "reads", "sortable": True},
                                            {"name": "action", "label": "View in IGV", "field": "action", "sortable": False}
                                        ]
                                        
                                        # Format rows for display
                                        rows = []
                                        fusion_region_map = {}  # Store region info for each row
                                        
                                        # Import re for position parsing
                                        import re
                                        
                                        for idx, row in clustered_data.iterrows():
                                            # Try to get start/end coordinates from the row if available
                                            # (e.g., from breakpoint validation)
                                            if all(col in row for col in ["gene1_start", "gene1_end", "gene2_start", "gene2_end"]):
                                                # Use actual start/end coordinates if available
                                                start1_raw = int(row["gene1_start"])
                                                end1_raw = int(row["gene1_end"])
                                                start2_raw = int(row["gene2_start"])
                                                end2_raw = int(row["gene2_end"])
                                            else:
                                                # Parse from position strings (format: "start-end")
                                                pos1_str = str(row.get("gene1_position", ""))
                                                pos2_str = str(row.get("gene2_position", ""))
                                                
                                                # Parse range format "start-end" or just a single number
                                                pos1_match = re.match(r'(\d+)[-–—](\d+)', pos1_str.replace(',', ''))
                                                pos2_match = re.match(r'(\d+)[-–—](\d+)', pos2_str.replace(',', ''))
                                                
                                                if pos1_match and pos2_match:
                                                    # Extract both start and end from the range
                                                    start1_raw = int(pos1_match.group(1))
                                                    end1_raw = int(pos1_match.group(2))
                                                    start2_raw = int(pos2_match.group(1))
                                                    end2_raw = int(pos2_match.group(2))
                                                else:
                                                    # Fallback: try to extract single coordinates
                                                    pos1_single = re.search(r'(\d+)', pos1_str.replace(',', ''))
                                                    pos2_single = re.search(r'(\d+)', pos2_str.replace(',', ''))
                                                    if pos1_single and pos2_single:
                                                        # Single coordinate - use it as both start and end
                                                        start1_raw = end1_raw = int(pos1_single.group(1))
                                                        start2_raw = end2_raw = int(pos2_single.group(1))
                                                    else:
                                                        # Skip this row if we can't parse coordinates
                                                        continue
                                            
                                            # Use the actual range (start to end), then add/subtract 10kb padding
                                            padding = 10000
                                            min1 = min(start1_raw, end1_raw)
                                            max1 = max(start1_raw, end1_raw)
                                            min2 = min(start2_raw, end2_raw)
                                            max2 = max(start2_raw, end2_raw)
                                            
                                            # Subtract 10kb from min and add 10kb to max
                                            start1 = max(1, min1 - padding)
                                            end1 = max1 + padding
                                            start2 = max(1, min2 - padding)
                                            end2 = max2 + padding
                                            
                                            chr1 = str(row.get("chr1", "Unknown"))
                                            chr2 = str(row.get("chr2", "Unknown"))
                                            
                                            # Format region as "chr1:start-end chr2:start-end"
                                            region = f"{chr1}:{start1}-{end1} {chr2}:{start2}-{end2}"
                                            
                                            # For display, show the breakpoint range (not the padded version)
                                            # Use a single representative coordinate or the range midpoint
                                            display_pos1 = f"{min1:,}-{max1:,}" if min1 != max1 else f"{min1:,}"
                                            display_pos2 = f"{min2:,}-{max2:,}" if min2 != max2 else f"{min2:,}"
                                            
                                            formatted_row = {
                                                "fusion_pair": row.get("fusion_pair", ""),
                                                "chr1": chr1,
                                                "pos1": display_pos1,
                                                "chr2": chr2,
                                                "pos2": display_pos2,
                                                "reads": int(row.get("reads", 0)),
                                                "action": "",
                                            }
                                            
                                            rows.append(formatted_row)
                                            fusion_region_map[len(rows) - 1] = region
                                        
                                        if rows:
                                            # Store fusion regions mapped by fusion pair for easy lookup
                                            fusion_regions_by_pair = {}
                                            for idx, row_data in enumerate(rows):
                                                fusion_regions_by_pair[row_data["fusion_pair"]] = fusion_region_map[idx]
                                            
                                            # Create JavaScript map of regions for IGV navigation
                                            import json
                                            js_regions_json = json.dumps(fusion_regions_by_pair)
                                            
                                            # Function to navigate IGV to a fusion region
                                            def navigate_to_fusion_region(fusion_pair: str):
                                                """Navigate IGV browser to the specified fusion pair region."""
                                                if fusion_pair in fusion_regions_by_pair:
                                                    region = fusion_regions_by_pair[fusion_pair]
                                                    # Escape region string for JavaScript
                                                    escaped_region = region.replace('"', '\\"').replace("'", "\\'")
                                                    js_navigate = f"""
                                                        (function() {{
                                                            try {{
                                                                if (window.lj_igv && window.lj_igv_browser_ready) {{
                                                                    console.log('[IGV] Navigating to fusion region: {escaped_region}');
                                                                    window.lj_igv.search('{escaped_region}');
                                                                }} else {{
                                                                    console.warn('[IGV] Browser not ready yet, will navigate when ready');
                                                                    setTimeout(function() {{
                                                                        if (window.lj_igv && window.lj_igv_browser_ready) {{
                                                                            window.lj_igv.search('{escaped_region}');
                                                                        }}
                                                                    }}, 500);
                                                                }}
                                                            }} catch (error) {{
                                                                console.error('[IGV] Navigation error:', error);
                                                            }}
                                                        }})();
                                                    """
                                                    ui.run_javascript(js_navigate, timeout=5.0)
                                            
                                            # Initialize fusion regions map in window BEFORE creating table
                                            # This ensures it's available when the slot template renders
                                            js_init_regions = f"""
                                                (function() {{
                                                    window.fusionRegionsMap = window.fusionRegionsMap || {{}};
                                                    Object.assign(window.fusionRegionsMap, {js_regions_json});
                                                    console.log('[Fusion] Initialized fusion regions map with', Object.keys(window.fusionRegionsMap).length, 'regions');
                                                }})();
                                            """
                                            ui.run_javascript(js_init_regions, timeout=5.0)
                                            
                                            # Create styled table
                                            table_container, fusion_table = styled_table(
                                                columns=columns,
                                                rows=rows,
                                                pagination=20,
                                                class_size="table-xs"
                                            )
                                            
                                            # Add clickable action button column using slot that emits events to Python
                                            try:
                                                # Use Vue event emission which is more reliable than window functions
                                                fusion_table.add_slot(
                                                    "body-cell-action",
                                                    """
<q-td key="action" :props="props">
  <q-btn 
    icon="visibility" 
    size="sm" 
    dense 
    flat 
    color="primary"
    @click="$parent.$emit('fusion-view-igv', props.row.fusion_pair)"
    title="View in IGV"
  />
</q-td>
"""
                                                )
                                                
                                                # Handle the event from the slot
                                                def on_fusion_view_igv(e):
                                                    """Handle fusion view IGV event from table button."""
                                                    try:
                                                        fusion_pair = e.args if isinstance(e.args, str) else getattr(e, 'args', None)
                                                        if fusion_pair:
                                                            logging.debug(f"[Fusion] Button clicked for: {fusion_pair}")
                                                            navigate_to_fusion_region(fusion_pair)
                                                    except Exception as ex:
                                                        logging.warning(f"Error handling fusion view IGV event: {ex}")
                                                
                                                fusion_table.on("fusion-view-igv", on_fusion_view_igv)
                                                logging.debug("Added action button column slot with event handler")
                                            except Exception as slot_ex:
                                                logging.warning(f"Could not add action column slot: {slot_ex}")
                                                import traceback
                                                logging.warning(traceback.format_exc())
                                                
                                                # Fallback: Use JavaScript with inline region lookup
                                                try:
                                                    # Create a safer inline handler that doesn't rely on window functions
                                                    js_inline_handler = f"""
                                                        (function() {{
                                                            // Store regions as a constant in the closure
                                                            const fusionRegions = {js_regions_json};
                                                            
                                                            // Create handler function immediately
                                                            window.handleFusionNav = function(fusionPair) {{
                                                                console.log('[Fusion] Navigation handler called for:', fusionPair);
                                                                const region = fusionRegions[fusionPair];
                                                                if (region) {{
                                                                    console.log('[Fusion] Navigating to:', region);
                                                                    function nav() {{
                                                                        if (window.lj_igv && window.lj_igv_browser_ready) {{
                                                                            window.lj_igv.search(region);
                                                                            return true;
                                                                        }}
                                                                        return false;
                                                                    }}
                                                                    if (!nav()) {{
                                                                        setTimeout(nav, 500);
                                                                        setTimeout(nav, 1500);
                                                                    }}
                                                                }}
                                                            }};
                                                            console.log('[Fusion] Created navigation handler');
                                                        }})();
                                                    """
                                                    ui.run_javascript(js_inline_handler, timeout=5.0)
                                                    
                                                    # Wait a moment for JS to execute, then add slot
                                                    ui.timer(0.1, lambda: None, once=True)
                                                    
                                                    fusion_table.add_slot(
                                                        "body-cell-action",
                                                        """
<q-td key="action" :props="props">
  <q-btn 
    icon="visibility" 
    size="sm" 
    dense 
    flat 
    color="primary"
    @click="window.handleFusionNav && window.handleFusionNav(props.row.fusion_pair)"
    title="View in IGV"
  />
</q-td>
"""
                                                    )
                                                except Exception as fallback_ex:
                                                    logging.warning(f"Fallback approach also failed: {fallback_ex}")
                                            
                                            # Use JavaScript to attach click handlers to rows - improved with event delegation
                                            js_attach_handlers = f"""
                                                (function() {{
                                                    let attached = false;
                                                    
                                                    function attachFusionHandlers() {{
                                                        if (attached) return;
                                                        
                                                        console.log('[Fusion] Attaching click handlers...');
                                                        
                                                        // Use event delegation on the document body for better reliability
                                                        function handleFusionRowClick(e) {{
                                                            // Check if click is on a fusion table row
                                                            let target = e.target;
                                                            let row = target.closest('tbody tr');
                                                            
                                                            if (!row) return;
                                                            
                                                            // Check if this row is in a fusion table
                                                            const table = row.closest('table');
                                                            if (!table) return;
                                                            
                                                            const headers = table.querySelectorAll('thead th');
                                                            let isFusionTable = false;
                                                            for (let i = 0; i < headers.length; i++) {{
                                                                if ((headers[i].textContent || '').toLowerCase().includes('fusion pair')) {{
                                                                    isFusionTable = true;
                                                                    break;
                                                                }}
                                                            }}
                                                            
                                                            if (!isFusionTable) return;
                                                            
                                                            // Skip if clicking on the action button (let button handle it)
                                                            if (target.closest('button') || target.closest('.q-btn')) {{
                                                                return;
                                                            }}
                                                            
                                                            e.preventDefault();
                                                            e.stopPropagation();
                                                            
                                                            console.log('[Fusion] Row clicked');
                                                            
                                                            // Get fusion pair from first cell
                                                            const cells = row.querySelectorAll('td');
                                                            if (cells.length > 0) {{
                                                                const fusionPair = cells[0].textContent.trim();
                                                                console.log('[Fusion] Fusion pair:', fusionPair);
                                                                
                                                                if (fusionPair && window.fusionRegionsMap && window.fusionRegionsMap[fusionPair]) {{
                                                                    const region = window.fusionRegionsMap[fusionPair];
                                                                    console.log('[Fusion] Navigating to region:', region);
                                                                    
                                                                    // Visual feedback
                                                                    row.style.backgroundColor = '#e3f2fd';
                                                                    setTimeout(function() {{
                                                                        row.style.backgroundColor = '';
                                                                    }}, 400);
                                                                    
                                                                    // Navigate IGV
                                                                    function navIGV() {{
                                                                        if (window.lj_igv && window.lj_igv_browser_ready) {{
                                                                            try {{
                                                                                window.lj_igv.search(region);
                                                                                console.log('[Fusion] IGV navigation successful');
                                                                                return true;
                                                                            }} catch(err) {{
                                                                                console.error('[Fusion] IGV error:', err);
                                                                                return false;
                                                                            }}
                                                                        }}
                                                                        return false;
                                                                    }}
                                                                    
                                                                    if (!navIGV()) {{
                                                                        setTimeout(function() {{ navIGV(); }}, 500);
                                                                        setTimeout(function() {{ navIGV(); }}, 1500);
                                                                    }}
                                                                }}
                                                            }}
                                                        }}
                                                        
                                                        // Attach event listener to document (event delegation)
                                                        document.addEventListener('click', handleFusionRowClick, true);
                                                        
                                                        // Also make rows visually clickable
                                                        const tables = document.querySelectorAll('table');
                                                        tables.forEach(function(table) {{
                                                            const headers = table.querySelectorAll('thead th');
                                                            let isFusionTable = false;
                                                            for (let i = 0; i < headers.length; i++) {{
                                                                if ((headers[i].textContent || '').toLowerCase().includes('fusion pair')) {{
                                                                    isFusionTable = true;
                                                                    break;
                                                                }}
                                                            }}
                                                            
                                                            if (isFusionTable) {{
                                                                const tbody = table.querySelector('tbody');
                                                                if (tbody) {{
                                                                    const rows = tbody.querySelectorAll('tr');
                                                                    rows.forEach(function(row) {{
                                                                        row.style.cursor = 'pointer';
                                                                    }});
                                                                    console.log('[Fusion] Made', rows.length, 'rows clickable');
                                                                }}
                                                            }}
                                                        }});
                                                        
                                                        attached = true;
                                                        console.log('[Fusion] Click handlers attached successfully');
                                                    }}
                                                    
                                                    // Try multiple times to ensure table is rendered
                                                    attachFusionHandlers();
                                                    setTimeout(attachFusionHandlers, 300);
                                                    setTimeout(attachFusionHandlers, 800);
                                                    setTimeout(attachFusionHandlers, 1500);
                                                    setTimeout(attachFusionHandlers, 2500);
                                                }})();
                                            """
                                            
                                            # Execute JavaScript after table is created
                                            ui.timer(0.5, lambda: ui.run_javascript(js_attach_handlers, timeout=10.0), once=True)
                                            ui.timer(2.0, lambda: ui.run_javascript(js_attach_handlers, timeout=10.0), once=True)
                                            
                                            # Add summary
                                            total_fusions = len(rows)
                                            total_reads = sum(r["reads"] for r in rows)
                                            ui.label(
                                                f"Total fusions: {total_fusions} | "
                                                f"Total supporting reads: {total_reads}"
                                            ).classes("classification-insight-foot")
                                        else:
                                            ui.label("No fusion pairs found").classes(
                                                "classification-insight-meta"
                                            )
                                else:
                                    with ui.element("div").classes(
                                        "classification-insight-shell w-full min-w-0"
                                    ):
                                        ui.label("Fusion pairs").classes(
                                            "classification-insight-heading text-headline-small"
                                        )
                                        ui.label(
                                            "No validated fusion pairs found in this run."
                                        ).classes("classification-insight-meta")
                            else:
                                with ui.element("div").classes(
                                    "classification-insight-shell w-full min-w-0"
                                ):
                                    ui.label("Fusion pairs").classes(
                                        "classification-insight-heading text-headline-small"
                                    )
                                    ui.label("Fusion data is empty.").classes(
                                        "classification-insight-meta"
                                    )
                        else:
                            with ui.element("div").classes(
                                "classification-insight-shell w-full min-w-0"
                            ):
                                ui.label("Fusion pairs").classes(
                                    "classification-insight-heading text-headline-small"
                                )
                                ui.label(
                                    "Fusion analysis data not available. Run fusion analysis first."
                                ).classes("classification-insight-meta")
                    
                    # Target Genes Table section
                    if sample_dir and sample_dir.exists():
                        target_coverage_file = sample_dir / "target_coverage.csv"
                        bed_coverage_file = sample_dir / "bed_coverage_main.csv"
                        
                        # Try target_coverage.csv first, fallback to bed_coverage_main.csv
                        coverage_file = None
                        if target_coverage_file.exists():
                            coverage_file = target_coverage_file
                        elif bed_coverage_file.exists():
                            coverage_file = bed_coverage_file
                        
                        if coverage_file:
                            # Load coverage data asynchronously to avoid blocking GUI
                            def load_coverage_data():
                                try:
                                    import pandas as pd
                                    return pd.read_csv(coverage_file)
                                except Exception as e:
                                    logging.error(f"Failed to load coverage file {coverage_file}: {e}")
                                    return None
                            
                            # Create placeholder table that will be updated when data loads
                            coverage_table_placeholder = ui.table(
                                columns=[
                                    {"name": "gene", "label": "Gene", "field": "gene"},
                                    {"name": "chrom", "label": "Chrom", "field": "chrom"},
                                    {"name": "startpos", "label": "Start", "field": "startpos"},
                                    {"name": "endpos", "label": "End", "field": "endpos"},
                                    {"name": "coverage", "label": "Coverage", "field": "coverage"},
                                ],
                                rows=[],
                            ).classes("w-full")
                            
                            # Load coverage data synchronously (fast enough for typical file sizes)
                            # Note: This runs during page creation, not during user interaction
                            try:
                                import pandas as pd
                                df = pd.read_csv(coverage_file)
                                
                                # Ensure we have the required columns
                                required_cols = ["chrom", "startpos", "endpos", "name"]
                                if all(col in df.columns for col in required_cols):
                                    # Calculate coverage if not present
                                    if "coverage" not in df.columns:
                                        if "length" in df.columns and "bases" in df.columns:
                                            df["coverage"] = df["bases"] / df["length"]
                                        elif "startpos" in df.columns and "endpos" in df.columns:
                                            df["length"] = df["endpos"] - df["startpos"] + 1
                                            if "bases" in df.columns:
                                                df["coverage"] = df["bases"] / df["length"]
                                            else:
                                                df["coverage"] = 0
                                    
                                    # Prepare table data
                                    table_data = []
                                    gene_regions_by_name = {}  # Store regions for navigation
                                    for _, row in df.iterrows():
                                        gene_name = str(row["name"])
                                        chrom = str(row["chrom"])
                                        startpos_raw = int(row["startpos"])
                                        endpos_raw = int(row["endpos"])
                                        
                                        # Store region info for navigation (with 10kb padding)
                                        padding = 10000
                                        startpos_nav = max(1, startpos_raw - padding)
                                        endpos_nav = endpos_raw + padding
                                        region = f"{chrom}:{startpos_nav}-{endpos_nav}"
                                        gene_regions_by_name[gene_name] = region
                                        
                                        table_data.append({
                                            "chrom": chrom,
                                            "startpos": f"{startpos_raw:,}",  # Format with commas for display
                                            "endpos": f"{endpos_raw:,}",  # Format with commas for display
                                            "name": gene_name,
                                            "coverage": float(row.get("coverage", 0)),
                                            "action": "",
                                        })
                                    
                                    if table_data:
                                        with ui.element("div").classes(
                                            "classification-insight-shell w-full min-w-0"
                                        ):
                                            ui.label("Target genes").classes(
                                                "classification-insight-heading text-headline-small"
                                            )
                                            ui.label(
                                                "Click a gene row to open the region in IGV."
                                            ).classes("classification-insight-meta w-full mb-2")
                                            
                                            from robin.gui.theme import styled_table
                                            
                                            columns = [
                                                {"name": "name", "label": "Gene Name", "field": "name", "sortable": True},
                                                {"name": "chrom", "label": "Chromosome", "field": "chrom", "sortable": True},
                                                {"name": "startpos", "label": "Start", "field": "startpos", "sortable": True},
                                                {"name": "endpos", "label": "End", "field": "endpos", "sortable": True},
                                                {"name": "coverage", "label": "Coverage (x)", "field": "coverage", "sortable": True},
                                                {"name": "action", "label": "View in IGV", "field": "action", "sortable": False}
                                            ]
                                            
                                            table_container, gene_table = styled_table(
                                                columns=columns,
                                                rows=table_data,
                                                pagination=25,
                                                class_size="table-xs"
                                            )
                                            
                                            # Add search functionality
                                            try:
                                                with gene_table.add_slot("top-right"):
                                                    with ui.input(placeholder="Search genes...").props("type=search").bind_value(
                                                        gene_table, "filter"
                                                    ).add_slot("append"):
                                                        ui.icon("search")
                                            except Exception:
                                                pass
                                            
                                            # Add colored coverage badges
                                            try:
                                                gene_table.add_slot(
                                                    "body-cell-coverage",
                                                    """
    <q-td key="coverage" :props="props">
    <q-badge :color="props.value >= 30 ? 'green' : props.value >= 20 ? 'blue' : props.value >= 10 ? 'orange' : 'red'">
        {{ Number(props.value).toFixed(2) }}
    </q-badge>
    </q-td>
    """,
                                                )
                                            except Exception:
                                                pass
                                            
                                            # Function to navigate IGV to a gene region
                                            import json
                                            js_gene_regions_json = json.dumps(gene_regions_by_name)
                                            
                                            def navigate_to_gene_region(gene_name: str):
                                                """Navigate IGV browser to the specified gene region."""
                                                if gene_name in gene_regions_by_name:
                                                    region = gene_regions_by_name[gene_name]
                                                    # Escape region string for JavaScript
                                                    escaped_region = region.replace('"', '\\"').replace("'", "\\'")
                                                    js_navigate = f"""
                                                        (function() {{
                                                            try {{
                                                                if (window.lj_igv && window.lj_igv_browser_ready) {{
                                                                    console.log('[IGV] Navigating to gene region: {escaped_region}');
                                                                    window.lj_igv.search('{escaped_region}');
                                                                }} else {{
                                                                    console.warn('[IGV] Browser not ready yet, will navigate when ready');
                                                                    setTimeout(function() {{
                                                                        if (window.lj_igv && window.lj_igv_browser_ready) {{
                                                                            window.lj_igv.search('{escaped_region}');
                                                                        }}
                                                                    }}, 500);
                                                                }}
                                                            }} catch (error) {{
                                                                console.error('[IGV] Navigation error:', error);
                                                            }}
                                                        }})();
                                                    """
                                                    ui.run_javascript(js_navigate, timeout=5.0)
                                            
                                            # Store regions in window object for JavaScript access
                                            js_init_gene_regions = f"""
                                                window.geneRegionsMap = window.geneRegionsMap || {{}};
                                                Object.assign(window.geneRegionsMap, {js_gene_regions_json});
                                                console.log('[Gene] Loaded', Object.keys(window.geneRegionsMap).length, 'gene regions');
                                            """
                                            ui.run_javascript(js_init_gene_regions, timeout=5.0)
                                            
                                            # Add clickable action button column using slot that emits events to Python
                                            try:
                                                # Use Vue event emission which is more reliable than window functions
                                                gene_table.add_slot(
                                                    "body-cell-action",
                                                    """
<q-td key="action" :props="props">
  <q-btn 
    icon="visibility" 
    size="sm" 
    dense 
    flat 
    color="primary"
    @click="$parent.$emit('gene-view-igv', props.row.name)"
    title="View in IGV"
  />
</q-td>
"""
                                                )
                                                
                                                # Handle the event from the slot
                                                def on_gene_view_igv(e):
                                                    """Handle gene view IGV event from table button."""
                                                    try:
                                                        gene_name = e.args if isinstance(e.args, str) else getattr(e, 'args', None)
                                                        if gene_name:
                                                            logging.debug(f"[Gene] Button clicked for: {gene_name}")
                                                            navigate_to_gene_region(gene_name)
                                                    except Exception as ex:
                                                        logging.warning(f"Error handling gene view IGV event: {ex}")
                                                
                                                gene_table.on("gene-view-igv", on_gene_view_igv)
                                                logging.debug("Added action button column slot with event handler")
                                            except Exception as slot_ex:
                                                logging.warning(f"Could not add action column slot: {slot_ex}")
                                                import traceback
                                                logging.warning(traceback.format_exc())
                                            
                                            # Also add row click handlers as a fallback
                                            try:
                                                # Create JavaScript to handle row clicks
                                                js_gene_table_handlers = f"""
                                                    (function() {{
                                                        // Store regions in window for JavaScript access
                                                        window.geneRegionsMap = window.geneRegionsMap || {{}};
                                                        Object.assign(window.geneRegionsMap, {js_gene_regions_json});
                                                        
                                                        function attachGeneTableHandlers() {{
                                                            // Find tables with "Gene Name" header
                                                            const tables = document.querySelectorAll('table');
                                                            
                                                            tables.forEach(function(table) {{
                                                                const headers = table.querySelectorAll('thead th');
                                                                let isGeneTable = false;
                                                                
                                                                for (let i = 0; i < headers.length; i++) {{
                                                                    const headerText = (headers[i].textContent || '').toLowerCase();
                                                                    if (headerText.includes('gene name') && headerText.includes('chromosome')) {{
                                                                        isGeneTable = true;
                                                                        break;
                                                                    }}
                                                                }}
                                                                
                                                                if (isGeneTable) {{
                                                                    const tbody = table.querySelector('tbody');
                                                                    if (tbody) {{
                                                                        const rows = tbody.querySelectorAll('tr');
                                                                        rows.forEach(function(row) {{
                                                                            // Skip if clicking on the action button (let button handle it)
                                                                            if (row.hasAttribute('data-gene-handled')) return;
                                                                            row.setAttribute('data-gene-handled', 'true');
                                                                            
                                                                            // Make row clickable (but button takes precedence)
                                                                            row.style.cursor = 'pointer';
                                                                            
                                                                            row.onclick = function(e) {{
                                                                                // Don't handle if clicking on button
                                                                                if (e.target.closest('button') || e.target.closest('.q-btn')) {{
                                                                                    return;
                                                                                }}
                                                                                
                                                                                e.preventDefault();
                                                                                e.stopPropagation();
                                                                                
                                                                                const cells = row.querySelectorAll('td');
                                                                                if (cells.length >= 5) {{
                                                                                    const geneName = cells[0].textContent.trim();
                                                                                    
                                                                                    if (geneName && window.geneRegionsMap && window.geneRegionsMap[geneName]) {{
                                                                                        const region = window.geneRegionsMap[geneName];
                                                                                        console.log('[Gene] Row clicked, navigating to region:', region);
                                                                                        
                                                                                        // Visual feedback
                                                                                        const origBg = row.style.backgroundColor;
                                                                                        row.style.backgroundColor = '#e3f2fd';
                                                                                        setTimeout(function() {{
                                                                                            row.style.backgroundColor = origBg;
                                                                                        }}, 400);
                                                                                        
                                                                                        // Navigate IGV
                                                                                        function navIGV() {{
                                                                                            if (window.lj_igv && window.lj_igv_browser_ready) {{
                                                                                                try {{
                                                                                                    window.lj_igv.search(region);
                                                                                                    console.log('[Gene] IGV navigation successful');
                                                                                                    return true;
                                                                                                }} catch(err) {{
                                                                                                    console.error('[Gene] IGV error:', err);
                                                                                                    return false;
                                                                                                }}
                                                                                            }}
                                                                                            return false;
                                                                                        }}
                                                                                        
                                                                                        if (!navIGV()) {{
                                                                                            setTimeout(function() {{ navIGV(); }}, 500);
                                                                                            setTimeout(function() {{ navIGV(); }}, 1500);
                                                                                        }}
                                                                                    }}
                                                                                }}
                                                                            }};
                                                                        }});
                                                                    }}
                                                                }}
                                                            }});
                                                        }}
                                                        
                                                        attachGeneTableHandlers();
                                                        setTimeout(attachGeneTableHandlers, 300);
                                                        setTimeout(attachGeneTableHandlers, 800);
                                                        setTimeout(attachGeneTableHandlers, 1500);
                                                    }})();
                                                """
                                                
                                                ui.timer(0.5, lambda: ui.run_javascript(js_gene_table_handlers, timeout=10.0), once=True)
                                                ui.timer(2.0, lambda: ui.run_javascript(js_gene_table_handlers, timeout=10.0), once=True)
                                            except Exception as e:
                                                logging.warning(f"Could not add gene table click handlers: {e}")
                                            
                                            # Add summary
                                            total_genes = len(table_data)
                                            ui.label(
                                                f"Total target genes: {total_genes}"
                                            ).classes("classification-insight-foot")
                            except Exception as e:
                                logging.warning(f"Could not load target gene table: {e}")

    def _create_workflow_monitor(self):
        """Create the main workflow monitoring page."""
        with theme.frame(
            "R.O.B.I.N - Workflow Monitor",
            smalltitle="Samples",
            batphone=False,
            center=self.center,
            setup_notifications=self._setup_notification_system,
        ):
            with ui.element("div").classes("w-full min-w-0").props(
                "id=workflow-monitor-page"
            ):
                with ui.column().classes(
                    "w-full gap-3 p-2 md:p-3 max-w-6xl mx-auto"
                ):
                    with ui.element("div").classes(
                        "classification-insight-shell w-full min-w-0"
                    ):
                        ui.label("Workflow monitor").classes(
                            "classification-insight-heading text-headline-small"
                        )
                        ui.label(
                            "Real-time workflow monitoring and control."
                        ).classes("classification-insight-foot")

                    # Workflow status overview
                    with ui.element("div").classes(
                        "classification-insight-card w-full min-w-0"
                    ):
                        with ui.column().classes(
                            "w-full min-w-0 gap-3 p-2 md:p-3"
                        ):
                            with ui.row().classes("items-center gap-2 min-w-0"):
                                ui.icon("monitor_heart").classes(
                                    "classification-insight-icon"
                                )
                                ui.label("Workflow status").classes(
                                    "classification-insight-model flex-1 min-w-0"
                                )

                            with ui.row().classes("w-full items-center gap-3 min-w-0"):
                                self.status_indicator = ui.icon(
                                    "fiber_manual_record"
                                ).classes(
                                    "text-2xl workflow-monitor-status-icon "
                                    "workflow-monitor-status-icon--running"
                                )
                                self.status_label = ui.label(
                                    "Workflow status: Running"
                                ).classes(
                                    "text-sm font-medium "
                                    "workflow-monitor-status-text--running"
                                )

                            with ui.row().classes(
                                "w-full gap-2 mt-2 flex-wrap"
                            ):
                                self.workflow_start_time = ui.label(
                                    "Started: —"
                                ).classes("text-sm workflow-monitor-meta")
                                self.workflow_duration = ui.label(
                                    "Duration: —"
                                ).classes("text-sm workflow-monitor-meta")

                            self.progress_bar = (
                                ui.linear_progress(0.0)
                                .classes("w-full mt-4")
                                .style("color: transparent")
                            )
                            self.progress_label = ui.label("0% complete").classes(
                                "text-sm text-center w-full workflow-monitor-meta"
                            )

                            ui.label("Run counts").classes(
                                "target-coverage-panel__meta-label mt-2 mb-1"
                            )
                            with ui.row().classes(
                                "w-full gap-4 mt-1 flex-wrap"
                            ):
                                with ui.row().classes("items-center gap-2"):
                                    ui.label("Completed").classes(
                                        "text-xs workflow-monitor-meta"
                                    )
                                    self.completed_count = ui.label("0").classes(
                                        "text-xs font-semibold workflow-monitor-num"
                                    )
                                with ui.row().classes("items-center gap-2"):
                                    ui.label("Failed").classes(
                                        "text-xs workflow-monitor-meta"
                                    )
                                    self.failed_count = ui.label("0").classes(
                                        "text-xs font-semibold workflow-monitor-num"
                                    )
                                with ui.row().classes("items-center gap-2"):
                                    ui.label("Total").classes(
                                        "text-xs workflow-monitor-meta"
                                    )
                                    self.total_count = ui.label("0").classes(
                                        "text-xs font-semibold workflow-monitor-num"
                                    )

                    # File processing progress (per run)
                    with ui.element("div").classes(
                        "classification-insight-card w-full min-w-0"
                    ):
                        with ui.column().classes(
                            "w-full min-w-0 gap-3 p-2 md:p-3"
                        ):
                            with ui.row().classes("items-center gap-2 min-w-0"):
                                ui.icon("folder_open").classes(
                                    "classification-insight-icon"
                                )
                                ui.label("File processing (per run)").classes(
                                    "classification-insight-model flex-1 min-w-0"
                                )

                            with ui.row().classes(
                                "w-full items-center gap-3 min-w-0"
                            ):
                                ui.label("Overall files").classes(
                                    "text-sm font-medium workflow-monitor-meta shrink-0"
                                )
                                self.overall_files_progress = (
                                    ui.linear_progress(0.0).classes(
                                        "flex-1 min-w-0"
                                    )
                                )
                                self.overall_files_label = ui.label(
                                    "0/0 files processed"
                                ).classes(
                                    "text-sm min-w-[120px] workflow-monitor-meta"
                                )

                            ui.label("Per-sample progress").classes(
                                "target-coverage-panel__meta-label mt-2 mb-1"
                            )
                            with ui.scroll_area().classes("w-full").style(
                                "max-height: 300px;"
                            ):
                                self.sample_files_progress_container = ui.column().classes(
                                    "w-full gap-2 p-1 min-w-0"
                                )

                    # Queue status
                    with ui.element("div").classes(
                        "classification-insight-card w-full min-w-0"
                    ):
                        with ui.column().classes(
                            "w-full min-w-0 gap-3 p-2 md:p-3"
                        ):
                            with ui.row().classes("items-center gap-2 min-w-0"):
                                ui.icon("layers").classes(
                                    "classification-insight-icon"
                                )
                                ui.label("Queue status").classes(
                                    "classification-insight-model flex-1 min-w-0"
                                )

                            with ui.row().classes(
                                "w-full gap-2 flex-wrap items-stretch"
                            ):
                                with ui.element("div").classes(
                                    "workflow-monitor-queue-tile flex-1 min-w-[10rem]"
                                ):
                                    ui.label("Preprocessing").classes(
                                        "classification-insight-foot"
                                    )
                                    self.preprocessing_status = ui.label("0/0").classes(
                                        "text-2xl font-bold workflow-monitor-queue-num--prep"
                                    )

                                with ui.element("div").classes(
                                    "workflow-monitor-queue-tile flex-1 min-w-[10rem]"
                                ):
                                    with ui.row().classes(
                                        "items-center gap-1 min-w-0"
                                    ):
                                        ui.icon("science").classes(
                                            "text-base workflow-monitor-queue-icon"
                                        )
                                        ui.label("Analysis").classes(
                                            "classification-insight-foot"
                                        )
                                    self.analysis_status = ui.label("0/0").classes(
                                        "text-2xl font-bold workflow-monitor-queue-num--analysis"
                                    )

                                with ui.element("div").classes(
                                    "workflow-monitor-queue-tile flex-1 min-w-[10rem]"
                                ):
                                    ui.label("Classification").classes(
                                        "classification-insight-foot"
                                    )
                                    self.classification_status = ui.label("0/0").classes(
                                        "text-2xl font-bold workflow-monitor-queue-num--class"
                                    )

                                with ui.element("div").classes(
                                    "workflow-monitor-queue-tile flex-1 min-w-[10rem]"
                                ):
                                    ui.label("Other").classes(
                                        "classification-insight-foot"
                                    )
                                    self.other_status = ui.label("0/0").classes(
                                        "text-2xl font-bold workflow-monitor-queue-num--other"
                                    )

                            try:
                                if self._last_queue_status:
                                    self._update_queue_status(self._last_queue_status)
                            except Exception:
                                pass

                    # Active jobs
                    with ui.element("div").classes(
                        "classification-insight-card w-full min-w-0"
                    ):
                        with ui.column().classes(
                            "w-full min-w-0 gap-3 p-2 md:p-3"
                        ):
                            with ui.row().classes("items-center gap-2 min-w-0"):
                                ui.icon("work").classes(
                                    "classification-insight-icon"
                                )
                                ui.label("Active jobs").classes(
                                    "classification-insight-model flex-1 min-w-0"
                                )

                            with ui.row().classes(
                                "items-center gap-2 mb-2 flex-wrap"
                            ):
                                self.active_jobs_search = ui.input("Search…").props(
                                    "outlined dense clearable"
                                ).classes("min-w-[12rem] flex-1")
                                self.active_jobs_type_filter = (
                                    ui.select(
                                        options=["All"],
                                        value="All",
                                        label="Type",
                                    )
                                    .props("outlined dense clearable")
                                    .classes("w-40")
                                )
                                self.active_jobs_worker_filter = (
                                    ui.select(
                                        options=["All"],
                                        value="All",
                                        label="Worker",
                                    )
                                    .props("outlined dense clearable")
                                    .classes("w-40")
                                )

                            from robin.gui.theme import styled_table

                            _jobs_container, self.active_jobs_table = styled_table(
                                columns=[
                                    {
                                        "name": "job_id",
                                        "label": "Job ID",
                                        "field": "job_id",
                                        "sortable": True,
                                    },
                                    {
                                        "name": "job_type",
                                        "label": "Type",
                                        "field": "job_type",
                                        "sortable": True,
                                    },
                                    {
                                        "name": "filepath",
                                        "label": "File",
                                        "field": "filepath",
                                        "sortable": True,
                                    },
                                    {
                                        "name": "worker",
                                        "label": "Worker",
                                        "field": "worker",
                                        "sortable": True,
                                    },
                                    {
                                        "name": "duration",
                                        "label": "Duration",
                                        "field": "duration",
                                        "sortable": True,
                                    },
                                    {
                                        "name": "progress",
                                        "label": "Progress",
                                        "field": "progress",
                                        "sortable": True,
                                    },
                                ],
                                rows=[],
                                pagination=10,
                                class_size="table-xs",
                            )
                            try:
                                self.active_jobs_table.props(
                                    'multi-sort rows-per-page-options="[10,20,50,0]"'
                                )
                                self.active_jobs_search.bind_value(
                                    self.active_jobs_table, "filter"
                                )
                            except Exception:
                                pass
                            try:

                                def _on_active_jobs_filter_change(_=None):
                                    self._apply_active_jobs_filters_and_update()

                                self.active_jobs_type_filter.on(
                                    "update:model-value",
                                    _on_active_jobs_filter_change,
                                )
                                self.active_jobs_worker_filter.on(
                                    "update:model-value",
                                    _on_active_jobs_filter_change,
                                )
                            except Exception:
                                pass

                            self.no_active_jobs_label = ui.label(
                                "No active jobs at the moment."
                            ).classes("classification-insight-foot mt-2")

                    # Live logs
                    with ui.element("div").classes(
                        "classification-insight-card w-full min-w-0"
                    ):
                        with ui.column().classes(
                            "w-full min-w-0 gap-3 p-2 md:p-3"
                        ):
                            with ui.row().classes(
                                "w-full items-center justify-between gap-2 "
                                "flex-wrap min-w-0"
                            ):
                                with ui.row().classes(
                                    "items-center gap-2 min-w-0"
                                ):
                                    ui.icon("article").classes(
                                        "classification-insight-icon"
                                    )
                                    ui.label("Live logs").classes(
                                        "classification-insight-model flex-1 min-w-0"
                                    )
                                with ui.row().classes("gap-2 shrink-0"):
                                    ui.button(
                                        "Clear",
                                        on_click=self._clear_logs,
                                        icon="clear_all",
                                    ).props("flat no-caps outline")
                                    ui.button(
                                        "Export",
                                        on_click=lambda: self._export_logs(),
                                        icon="download",
                                    ).props("color=primary no-caps outline")

                            self.log_area = (
                                ui.textarea(
                                    "Workflow logs will appear here…"
                                )
                                .classes(
                                    "w-full h-40 workflow-monitor-log-area"
                                )
                                .props("outlined readonly dense")
                            )

                    # Configuration
                    with ui.element("div").classes(
                        "classification-insight-card w-full min-w-0"
                    ):
                        with ui.column().classes(
                            "w-full min-w-0 gap-3 p-2 md:p-3"
                        ):
                            with ui.row().classes("items-center gap-2 min-w-0"):
                                ui.icon("tune").classes(
                                    "classification-insight-icon"
                                )
                                ui.label("Workflow configuration").classes(
                                    "classification-insight-model flex-1 min-w-0"
                                )

                            with ui.grid(columns=2).classes(
                                "w-full gap-3 min-w-0"
                            ):
                                with ui.column().classes("min-w-0 gap-1"):
                                    ui.label("Monitored directory").classes(
                                        "target-coverage-panel__meta-label"
                                    )
                                    ui.label(
                                        self.monitored_directory or "Not specified"
                                    ).classes(
                                        "text-sm break-all font-mono "
                                        "workflow-monitor-config-value"
                                    )

                                    ui.label("Workflow steps").classes(
                                        "target-coverage-panel__meta-label mt-2"
                                    )
                                    ui.label(
                                        ", ".join(self.workflow_steps)
                                        if self.workflow_steps
                                        else "Not specified"
                                    ).classes(
                                        "text-sm workflow-monitor-config-value"
                                    )

                                with ui.column().classes("min-w-0 gap-1"):
                                    ui.label("Log level").classes(
                                        "target-coverage-panel__meta-label"
                                    )
                                    ui.label("—").classes(
                                        "text-sm workflow-monitor-config-value"
                                    )

                                    ui.label("Analysis workers").classes(
                                        "target-coverage-panel__meta-label mt-2"
                                    )
                                    ui.label("—").classes(
                                        "text-sm workflow-monitor-config-value"
                                    )

                    # Error summary
                    with ui.element("div").classes(
                        "classification-insight-card w-full min-w-0"
                    ):
                        with ui.column().classes(
                            "w-full min-w-0 gap-3 p-2 md:p-3"
                        ):
                            with ui.row().classes("items-center gap-2 min-w-0"):
                                ui.icon("error_outline").classes(
                                    "classification-insight-icon"
                                )
                                ui.label("Errors & troubleshooting").classes(
                                    "classification-insight-model flex-1 min-w-0"
                                )

                            with ui.row().classes(
                                "w-full justify-between gap-2 flex-wrap"
                            ):
                                with ui.column().classes(
                                    "text-center min-w-[5rem]"
                                ):
                                    self.preprocessing_errors = ui.label("0").classes(
                                        "text-xl font-bold workflow-monitor-error-num"
                                    )
                                    ui.label("Preprocessing").classes(
                                        "text-xs workflow-monitor-meta"
                                    )
                                with ui.column().classes(
                                    "text-center min-w-[5rem]"
                                ):
                                    self.analysis_errors = ui.label("0").classes(
                                        "text-xl font-bold workflow-monitor-error-num"
                                    )
                                    ui.label("Analysis").classes(
                                        "text-xs workflow-monitor-meta"
                                    )
                                with ui.column().classes(
                                    "text-center min-w-[5rem]"
                                ):
                                    self.classification_errors = ui.label("0").classes(
                                        "text-xl font-bold workflow-monitor-error-num"
                                    )
                                    ui.label("Classification").classes(
                                        "text-xs workflow-monitor-meta"
                                    )

                            ui.separator().classes("mgmt-detail-separator")
                            ui.label("Recent errors").classes(
                                "target-coverage-panel__meta-label mt-1 mb-1"
                            )
                            self.error_summary_label = ui.label(
                                "No errors detected."
                            ).classes(
                                "text-xs workflow-monitor-error-summary"
                            )

                    ui.label(
                        "R.O.B.I.N workflow monitor — session active"
                    ).classes(
                        "classification-insight-foot text-center w-full py-2"
                    )

                # Signal that GUI is ready to receive updates
                self.gui_ready.set()
                logging.info("[GUI] UI created and ready to receive updates")
                # Start a periodic UI-thread drain of the update queue
                app.timer(0.3, self._drain_updates_on_ui, active=True)
                # Start duration refresher
                app.timer(1.0, lambda: self._refresh_duration(), active=True)

    def _refresh_duration(self):
        if self._is_running and self._start_time and hasattr(self, "workflow_duration"):
            try:
                elapsed_seconds = int(time.time() - self._start_time)
                self.workflow_duration.set_text(
                    f"Duration: {self._format_duration(elapsed_seconds)}"
                )
            except Exception:
                pass

    def _export_logs(self):
        """Export logs to a file."""
        try:
            # Simple log export functionality
            log_content = "".join(self._log_buffer)
            if log_content:
                # Create a simple download
                ui.download(log_content, "workflow_logs.txt")
            else:
                ui.notify("No logs to export", type="warning")
        except Exception as e:
            ui.notify(f"Export failed: {e}", type="error")

    def _clear_logs(self):
        try:
            self._log_buffer.clear()
            self.log_area.set_value("")
        except Exception:
            pass

    def _locate_reference_for_sample(self, sample_dir: Path) -> Optional[str]:
        """Locate a reference genome for SNP analysis, preferring CLI-specified value."""

        # 1. Prefer reference supplied on the command line / workflow runner
        try:
            runner = getattr(self, "workflow_runner", None)
            runner_reference = None

            if runner is not None:
                runner_reference = getattr(runner, "reference", None)

                # Some wrappers may store the reference on an inner object
                if not runner_reference and hasattr(runner, "manager"):
                    runner_reference = getattr(runner.manager, "reference", None)

            if runner_reference:
                ref_path = Path(str(runner_reference))
                if ref_path.exists():
                    logging.info(
                        "Using reference genome specified via command line: %s",
                        ref_path,
                    )
                    return str(ref_path)
                logging.warning(
                    "Reference genome provided via command line does not exist: %s",
                    ref_path,
                )
        except Exception as exc:
            logging.debug(f"Error checking workflow runner for reference genome: {exc}")

        # 2. Fallback to legacy heuristics (sample directory, monitored directory, env var)
        candidates: List[Path] = []

        try:
            candidates.extend([
                sample_dir / "reference.fasta",
                sample_dir / "reference.fa",
            ])

            base_dir = Path(self.monitored_directory) if self.monitored_directory else sample_dir.parent
            if base_dir:
                candidates.extend([
                    base_dir / "reference.fasta",
                    base_dir / "reference.fa",
                ])

            env_reference = os.environ.get("robin_REFERENCE")
            if env_reference:
                candidates.append(Path(env_reference))
        except Exception as exc:
            logging.debug(f"Error assembling fallback reference candidates for {sample_dir}: {exc}")

        for candidate in candidates:
            try:
                if candidate and candidate.exists():
                    return str(candidate)
            except Exception:
                continue

        return None

    def _sample_has_snp_calling_outputs(self, sample_id: str) -> bool:
        """True when Clair3 / SNP pipeline outputs exist (same check as coverage UI)."""
        if not self.monitored_directory:
            return False
        clair = Path(self.monitored_directory) / sample_id / "clair3"
        snp_vcf = clair / "snpsift_output.vcf"
        indel_vcf = clair / "snpsift_indel_output.vcf"
        return snp_vcf.is_file() and indel_vcf.is_file()

    def _sample_snp_prerequisites_met(self, sample_dir: Path) -> tuple[bool, str]:
        """target.bam + non-empty targets_exceeding_threshold.bed (matches per-sample SNP UI)."""
        target_bam = sample_dir / "target.bam"
        targets_bed = sample_dir / "targets_exceeding_threshold.bed"
        if not target_bam.is_file():
            return False, "missing target.bam"
        if not targets_bed.is_file():
            return False, "missing targets_exceeding_threshold.bed"
        try:
            if not targets_bed.read_text(encoding="utf-8", errors="replace").strip():
                return False, "empty targets BED"
        except OSError as e:
            return False, f"cannot read targets BED: {e}"
        return True, ""

    def _list_samples_needing_snp_calling(self) -> List[str]:
        """Samples under work dir that are ready for SNP but lack Clair3 outputs."""
        logging.info(
            "[missing_snp_scan] thread=%s",
            threading.current_thread().name,
        )
        if not self.monitored_directory:
            return []
        base = Path(self.monitored_directory)
        if not base.is_dir():
            return []
        need: List[str] = []
        try:
            for d in sorted(base.iterdir()):
                if not d.is_dir():
                    continue
                sid = d.name
                if self._sample_has_snp_calling_outputs(sid):
                    continue
                ok, _ = self._sample_snp_prerequisites_met(d)
                if ok:
                    need.append(sid)
        except Exception as e:
            logging.debug(f"List samples needing SNP: {e}")
        return need

    def _is_mnpflex_enabled_for_gui(self) -> bool:
        """True when MNP-Flex credentials exist in the server environment."""
        username = os.getenv("MNPFLEX_USERNAME") or os.getenv("EPIGNOSTIX_USERNAME")
        password = os.getenv("MNPFLEX_PASSWORD") or os.getenv("EPIGNOSTIX_PASSWORD")
        return bool(username and password)

    def _mnpflex_results_dir_for_sample(
        self, sample_dir: Path, sample_id: str
    ) -> Optional[Path]:
        """Return directory containing `bundle_summary.json` (if results exist)."""
        try:
            candidates = [
                sample_dir / f"mnpflex_results_{sample_id}",
                sample_dir / "mnpflex_results",
                sample_dir,
            ]
            for candidate in candidates:
                if (candidate / "bundle_summary.json").exists():
                    return candidate
            for candidate in sorted(sample_dir.glob("mnpflex_results_*")):
                if (candidate / "bundle_summary.json").exists():
                    return candidate
        except Exception:
            return None
        return None

    def _mnpflex_parquet_path_for_sample(
        self, sample_dir: Path, sample_id: str
    ) -> Optional[Path]:
        """Preferred parquet is `<sample_id>.parquet`, else any `*.parquet`."""
        preferred = sample_dir / f"{sample_id}.parquet"
        if preferred.exists():
            return preferred
        try:
            matches = list(sample_dir.glob("*.parquet"))
            return matches[0] if matches else None
        except Exception:
            return None

    def _list_samples_needing_mnpflex_for_visible(
        self,
        rows_override: Optional[List[Dict[str, Any]]] = None,
    ) -> List[str]:
        """Eligible samples in the current table view that lack MNP-Flex results.

        `rows_override` allows the caller to pass a UI-thread-copied snapshot of rows
        (safer than touching NiceGUI objects from a background thread).
        """
        if not self.monitored_directory:
            return []
        if not self._is_mnpflex_enabled_for_gui():
            return []

        logging.info(
            "[missing_mnpflex_scan] thread=%s rows_override=%s",
            threading.current_thread().name,
            "provided" if rows_override is not None else "none",
        )

        if rows_override is not None:
            rows = rows_override
        else:
            try:
                rows = getattr(self.samples_table, "rows", None) or []
            except Exception:
                rows = []
            if not rows:
                rows = getattr(self, "_last_samples_rows", []) or []

        need: List[str] = []
        try:
            for r in rows:
                sid = str(r.get("sample_id") or "").strip()
                if not sid or sid == "unknown":
                    continue
                sample_dir = Path(self.monitored_directory) / sid
                if not sample_dir.exists():
                    continue

                if (
                    self._mnpflex_results_dir_for_sample(sample_dir, sid) is not None
                ):
                    continue

                try:
                    total_jobs = int(r.get("total_jobs") or 0)
                    completed_jobs = int(r.get("completed_jobs") or 0)
                    active_jobs = int(r.get("active_jobs") or 0)
                    pending_jobs = int(r.get("pending_jobs") or 0)
                except Exception:
                    continue

                # Only run once the sample workflow is fully complete/inactive.
                if total_jobs <= 0:
                    continue
                if not (
                    active_jobs == 0
                    and pending_jobs == 0
                    and completed_jobs >= total_jobs
                ):
                    continue

                if self._mnpflex_parquet_path_for_sample(sample_dir, sid) is None:
                    continue

                need.append(sid)
        except Exception as e:
            logging.debug(f"List samples needing MNP-Flex: {e}")

        return sorted(set(need))

    def _zip_paths_for_bulk_download(self, paths: List[str]) -> Optional[str]:
        """
        Safari blocks multiple programmatic downloads; package multiple files into one ZIP.
        Single file: return path unchanged.
        """
        valid = [p for p in paths if p and os.path.isfile(p)]
        if not valid:
            return None
        if len(valid) == 1:
            return valid[0]
        try:
            # Deterministic ZIP name: includes current date/time (no random mkstemp id).
            # If multiple exports happen in the same second, append a small counter.
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            base_name = f"robin_reports_{timestamp}"
            zpath = os.path.join(tempfile.gettempdir(), f"{base_name}.zip")
            if os.path.exists(zpath):
                i = 2
                while os.path.exists(zpath):
                    zpath = os.path.join(
                        tempfile.gettempdir(), f"{base_name}_{i}.zip"
                    )
                    i += 1
            used_names: Set[str] = set()
            with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as zf:
                for fp in valid:
                    base_name = Path(fp).name

                    # The export file names already include the sample id (e.g. "<sid>_run_report.pdf").
                    # Avoid prefixing with the parent folder name again, which caused duplicates like:
                    # "<sid>_<sid>_run_report.pdf".
                    arcname = base_name
                    if arcname in used_names:
                        # Collision: preserve both by adding a suffix.
                        i = 2
                        stem, ext = os.path.splitext(base_name)
                        while f"{stem}.{i}{ext}" in used_names:
                            i += 1
                        arcname = f"{stem}.{i}{ext}"
                    used_names.add(arcname)
                    zf.write(fp, arcname=arcname)
            return zpath
        except Exception as e:
            logging.error(f"Could not zip export files: {e}")
            return None

    def _wait_for_snp_outputs_or_timeout(
        self,
        sample_id: str,
        *,
        poll_s: float = 5.0,
        max_wait_s: float = 86400.0,
    ) -> bool:
        """Poll disk for SNP outputs after job submission (Ray submit returns before job finishes)."""
        deadline = time.time() + max_wait_s
        while time.time() < deadline:
            if self._sample_has_snp_calling_outputs(sample_id):
                return True
            time.sleep(poll_s)
        return False

    def _set_pipeline_status(
        self,
        sample_id: str,
        *,
        phase: str,
        progress: float,
        detail: str = "",
        level: str = "info",
    ) -> None:
        """Store lightweight per-sample finalize/SNP progress for GUI table display."""
        try:
            self._sample_pipeline_status[sample_id] = {
                "phase": phase,
                "progress": max(0.0, min(1.0, float(progress))),
                "detail": detail,
                "level": level,
                "updated_at": time.time(),
            }
        except Exception:
            pass

    def _wait_for_target_bam_or_timeout(
        self,
        sample_id: str,
        *,
        poll_s: float = 2.0,
        max_wait_s: float = 3600.0,
    ) -> bool:
        """Poll disk until `target.bam` and index exist for sample."""
        if not self.monitored_directory:
            return False
        sample_dir = Path(self.monitored_directory) / sample_id
        target_bam = sample_dir / "target.bam"
        target_bai = sample_dir / "target.bam.bai"
        deadline = time.time() + max_wait_s
        while time.time() < deadline:
            if target_bam.exists() and target_bai.exists():
                return True
            time.sleep(poll_s)
        return False

    def _resolve_submission_result(self, maybe_result: Any, operation: str) -> bool:
        """Normalize workflow submission return (sync bool or async coroutine)."""
        try:
            logging.debug(
                "workflow submit: %s return_type=%s",
                operation,
                type(maybe_result).__name__,
            )
            if asyncio.iscoroutine(maybe_result):
                logging.debug(
                    "workflow submit: %s running asyncio.run on coroutine", operation
                )
                result = bool(asyncio.run(maybe_result))
                logging.debug(
                    "workflow submit: %s coroutine result=%s", operation, result
                )
                return result
            result = bool(maybe_result)
            logging.debug("workflow submit: %s sync result=%s", operation, result)
            return result
        except Exception as e:
            logging.debug("workflow submit: %s exception: %s", operation, e)
            logging.error(f"{operation} failed: {e}", exc_info=True)
            return False

    def _bulk_snp_sequential_run(self, sample_ids: List[str]) -> None:
        """Run SNP analysis one sample at a time (sequential, like manual triggers)."""
        try:
            logging.info(
                "[bulk_snp] worker start thread=%s n=%d",
                threading.current_thread().name,
                len(sample_ids or []),
            )
            from robin.analysis.target_analysis import is_docker_available_for_snp_analysis

            docker_ok, docker_error = is_docker_available_for_snp_analysis()
            if not docker_ok:
                logging.warning(
                    "Bulk SNP skipped: Docker not available. %s", docker_error or ""
                )
                self.send_update(
                    UpdateType.WARNING_NOTIFICATION,
                    {
                        "title": "SNP batch skipped",
                        "message": f"Docker is required for SNP calling. {docker_error or ''}",
                        "level": "warning",
                    },
                    priority=6,
                )
                return

            workflow_runner = getattr(self, "workflow_runner", None)
            if not workflow_runner or not hasattr(
                workflow_runner, "submit_snp_analysis_job"
            ):
                logging.warning("Bulk SNP: no workflow_runner.submit_snp_analysis_job")
                self.send_update(
                    UpdateType.WARNING_NOTIFICATION,
                    {
                        "title": "SNP batch skipped",
                        "message": "Workflow runner is not available for SNP jobs.",
                        "level": "warning",
                    },
                    priority=6,
                )
                return

            total = len(sample_ids)
            for idx, sid in enumerate(sample_ids, start=1):
                self._set_pipeline_status(
                    sid,
                    phase="SNP preparation",
                    progress=0.7,
                    detail=f"Bulk queue item {idx}/{total}",
                )
                sample_dir = Path(self.monitored_directory) / sid
                ok, reason = self._sample_snp_prerequisites_met(sample_dir)
                if not ok:
                    self._set_pipeline_status(
                        sid,
                        phase="SNP skipped",
                        progress=0.7,
                        detail=str(reason)[:120],
                        level="warning",
                    )
                    logging.info("Bulk SNP skip %s: %s", sid, reason)
                    continue
                if self._sample_has_snp_calling_outputs(sid):
                    self._set_pipeline_status(
                        sid,
                        phase="SNP complete",
                        progress=1.0,
                        detail="Outputs already present",
                    )
                    continue

                reference_path = self._locate_reference_for_sample(sample_dir)
                if not reference_path:
                    self._set_pipeline_status(
                        sid,
                        phase="SNP skipped",
                        progress=0.7,
                        detail="No reference genome configured",
                        level="warning",
                    )
                    logging.warning("Bulk SNP skip %s: no reference genome", sid)
                    continue

                with self._snp_analysis_lock:
                    self._snp_analysis_running_sample = sid
                try:
                    logging.info(
                        "Bulk SNP (%d/%d): submitting %s", idx, total, sid
                    )
                    submitted = self._resolve_submission_result(
                        workflow_runner.submit_snp_analysis_job(
                            sample_dir=str(sample_dir),
                            sample_id=sid,
                            reference=reference_path,
                            threads=4,
                            force_regenerate=False,
                        ),
                        operation=f"Bulk SNP submission for {sid}",
                    )
                    if not submitted:
                        self._set_pipeline_status(
                            sid,
                            phase="SNP submit failed",
                            progress=0.8,
                            detail="Submission failed",
                            level="warning",
                        )
                        logging.warning("Bulk SNP: submit failed for %s", sid)
                        continue
                    self._set_pipeline_status(
                        sid,
                        phase="SNP queued",
                        progress=0.85,
                        detail="Queued in slow worker",
                    )
                    finished = self._wait_for_snp_outputs_or_timeout(
                        sid, poll_s=5.0, max_wait_s=86400.0
                    )
                    if finished:
                        self._set_pipeline_status(
                            sid,
                            phase="SNP complete",
                            progress=1.0,
                            detail="SNP outputs detected",
                        )
                        logging.info("Bulk SNP: completed %s", sid)
                    else:
                        self._set_pipeline_status(
                            sid,
                            phase="SNP still running",
                            progress=0.9,
                            detail="Timed out waiting for outputs",
                            level="warning",
                        )
                        logging.warning(
                            "Bulk SNP: timeout waiting for outputs for %s — continuing",
                            sid,
                        )
                finally:
                    with self._snp_analysis_lock:
                        self._snp_analysis_running_sample = None

            self.send_update(
                UpdateType.WARNING_NOTIFICATION,
                {
                    "title": "SNP batch finished",
                    "message": f"Processed queue of {total} sample(s). Check logs for any skips or timeouts.",
                    "level": "info",
                },
                priority=4,
            )
        except Exception as e:
            logging.error(f"Bulk SNP sequential run failed: {e}", exc_info=True)
            self.send_update(
                UpdateType.WARNING_NOTIFICATION,
                {
                    "title": "SNP batch error",
                    "message": str(e)[:500],
                    "level": "negative",
                },
                priority=6,
            )
        finally:
            self._bulk_snp_worker_running = False

    def _bulk_mnpflex_sequential_run(self, sample_ids: List[str]) -> None:
        """Run MNP-Flex sequentially for samples missing results."""
        try:
            logging.info(
                "[bulk_mnpflex] worker start thread=%s n=%d",
                threading.current_thread().name,
                len(sample_ids or []),
            )
            if not self._is_mnpflex_enabled_for_gui():
                self.send_update(
                    UpdateType.WARNING_NOTIFICATION,
                    {
                        "title": "MNP-Flex batch skipped",
                        "message": "Missing MNP-Flex credentials in the server environment.",
                        "level": "warning",
                    },
                    priority=6,
                )
                return

            username = os.getenv("MNPFLEX_USERNAME") or os.getenv("EPIGNOSTIX_USERNAME")
            password = os.getenv("MNPFLEX_PASSWORD") or os.getenv("EPIGNOSTIX_PASSWORD")
            base_url = os.getenv("MNPFLEX_BASE_URL", "https://app.epignostix.com")
            workflow_id_env = os.getenv("MNPFLEX_WORKFLOW_ID", "18")
            try:
                workflow_id = int(workflow_id_env)
            except ValueError:
                workflow_id = 18
                logging.warning(
                    "MNPFLEX_WORKFLOW_ID invalid (%s); defaulting to %s",
                    workflow_id_env,
                    workflow_id,
                )
            client_id = os.getenv("MNPFLEX_CLIENT_ID", "ROBIN")
            client_secret = os.getenv("MNPFLEX_CLIENT_SECRET", "SECRET")
            scope = os.getenv("MNPFLEX_SCOPE", "")

            from robin import resources as robin_resources
            from robin.analysis.utilities.matkit import (
                reconstruct_full_bedmethyl_for_mnpflex,
            )
            from robin.analysis.utilities.mnp_flex import APIClient as MnpFlexApiClient
            from robin.utils.mnpflex_client_standalone import MNPFlexClient

            reference_bed = os.path.join(
                os.path.dirname(os.path.abspath(robin_resources.__file__)),
                "mnp_flex_sample_clean.bed",
            )

            # Reuse streaming API client across samples.
            api_client = MnpFlexApiClient(base_url="https://mnp-flex.org", verify_ssl=False)

            total = len(sample_ids)
            processed = 0
            skipped = 0
            errors = 0

            for idx, sid in enumerate(sample_ids, start=1):
                sample_dir = (
                    Path(self.monitored_directory) / sid
                    if self.monitored_directory
                    else None
                )
                if not sample_dir or not sample_dir.exists():
                    skipped += 1
                    logging.warning("Bulk MNP-Flex skip %s: missing sample dir", sid)
                    continue

                if self._mnpflex_results_dir_for_sample(sample_dir, sid) is not None:
                    skipped += 1
                    continue

                parquet_path = self._mnpflex_parquet_path_for_sample(sample_dir, sid)
                if parquet_path is None:
                    skipped += 1
                    logging.warning("Bulk MNP-Flex skip %s: no parquet file", sid)
                    continue

                try:
                    logging.info(
                        "Bulk MNP-Flex (%d/%d): running %s", idx, total, sid
                    )

                    output_dir = sample_dir / f"mnpflex_results_{sid}"
                    output_dir.mkdir(parents=True, exist_ok=True)

                    # Build full bedmethyl BED from the parquet, then create the MNP-Flex subset BED.
                    bed_path = sample_dir / f"{sid}.mnpflex.bed"
                    bed_df = reconstruct_full_bedmethyl_for_mnpflex(str(parquet_path))
                    bed_df.to_csv(bed_path, sep="\t", index=False, header=False)

                    subset_path = sample_dir / f"{sid}.MNPFlex.subset.bed"
                    api_client.process_streaming(
                        reference_bed, str(bed_path), str(subset_path)
                    )

                    client = MNPFlexClient(
                        base_url=base_url,
                        username=username,
                        password=password,
                        verify_ssl=False,
                        client_id=client_id,
                        client_secret=client_secret,
                        scope=scope,
                    )
                    client.authenticate(
                        username=username,
                        password=password,
                        client_id=client_id,
                        client_secret=client_secret,
                    )
                    logging.info(
                        "Bulk MNP-Flex (%d/%d): submitting sample=%s workflow_id=%s subset_bed=%s output_dir=%s",
                        idx,
                        total,
                        sid,
                        workflow_id,
                        str(subset_path),
                        str(output_dir),
                    )
                    client.upload_retrieve_cleanup(
                        bed_file_path=str(subset_path),
                        sample_identifier=sid,
                        workflow_id=workflow_id,
                        output_dir=str(output_dir),
                    )
                    processed += 1
                except Exception as exc:
                    errors += 1
                    logging.error(
                        "Bulk MNP-Flex error for %s: %s", sid, exc, exc_info=True
                    )

            self.send_update(
                UpdateType.WARNING_NOTIFICATION,
                {
                    "title": "MNP-Flex batch finished",
                    "message": (
                        f"Processed {processed} sample(s), skipped {skipped}, errors {errors}. "
                        "Check logs for any failures."
                    ),
                    "level": "info",
                },
                priority=4,
            )
        except Exception as e:
            logging.error("Bulk MNP-Flex sequential run failed: %s", e, exc_info=True)
            self.send_update(
                UpdateType.WARNING_NOTIFICATION,
                {
                    "title": "MNP-Flex batch error",
                    "message": str(e)[:500],
                    "level": "negative",
                },
                priority=6,
            )
        finally:
            self._bulk_mnpflex_worker_running = False

    def _samples_select_all_visible_for_export(self) -> None:
        """Select every sample_id in the current (filtered) table rows for report export."""
        try:
            rows = getattr(self.samples_table, "rows", None) or []
            visible_ids = {str(r.get("sample_id")) for r in rows if r.get("sample_id")}
            self._selected_sample_ids = set(visible_ids)
            for r in self._last_samples_rows or []:
                sid = r.get("sample_id")
                if sid:
                    r["export"] = str(sid) in self._selected_sample_ids
            if self._selected_sample_ids:
                self.export_reports_button.enable()
            else:
                self.export_reports_button.disable()
            self._apply_samples_table_filters()
            logging.info(
                "[samples_overview] Select all (toolbar): %d sample(s) "
                "(visible filtered rows=%d)",
                len(self._selected_sample_ids),
                len(rows),
            )
        except Exception as e:
            logging.warning(
                "[samples_overview] Select all for export failed: %s", e, exc_info=True
            )

    def _samples_clear_export_selection(self) -> None:
        try:
            self._selected_sample_ids.clear()
            for r in self._last_samples_rows or []:
                r["export"] = False
            self.export_reports_button.disable()
            self._apply_samples_table_filters()
            logging.info("[samples_overview] Clear selection (toolbar)")
        except Exception as e:
            logging.warning(
                "[samples_overview] Clear export selection failed: %s",
                e,
                exc_info=True,
            )

    def _unlink_quiet(self, path: str) -> None:
        try:
            if path and os.path.isfile(path):
                os.remove(path)
        except OSError:
            pass

    def _trigger_target_bam_finalization(self, sample_id: str, *, trigger_snp: bool = False) -> None:
        """Trigger target.bam finalization for a sample. Optionally start SNP analysis afterwards."""
        logging.debug(
            "target_bam_finalize: start sample_id=%s trigger_snp=%s", sample_id, trigger_snp
        )
        self._set_pipeline_status(
            sample_id,
            phase="Finalize requested",
            progress=0.05,
            detail="Preparing finalization job",
        )
        if not self.monitored_directory:
            logging.debug("target_bam_finalize: no monitored_directory")
            self._set_pipeline_status(
                sample_id,
                phase="Finalize failed",
                progress=0.0,
                detail="No monitored directory available",
                level="negative",
            )
            return

        if self._is_target_bam_finalize_redundant(sample_id):
            self._finalized_samples.add(sample_id)
            logging.debug(
                "target_bam_finalize: redundant (target.bam already) sample_id=%s", sample_id
            )
            self._set_pipeline_status(
                sample_id,
                phase="Finalize already complete",
                progress=0.6 if trigger_snp else 1.0,
                detail="target.bam already present",
            )
            logging.info(
                f"Sample {sample_id} already has finalized target.bam on disk; "
                "skipping target.bam finalization"
            )
            if not trigger_snp:
                return
            logging.debug(
                "target_bam_finalize: redundant but trigger_snp=True, SNP path sample_id=%s",
                sample_id,
            )

        try:
            from robin.analysis.target_analysis import finalize_accumulation_for_sample
            import csv
            
            # Read target_panel from master.csv
            sample_dir = Path(self.monitored_directory) / sample_id
            master_csv = sample_dir / "master.csv"
            target_panel = "rCNS2"  # Default
            reference_path_for_finalize = self._locate_reference_for_sample(sample_dir)
            logging.debug(
                "target_bam_finalize: reference sample_id=%s path=%s",
                sample_id,
                reference_path_for_finalize,
            )

            if master_csv.exists():
                try:
                    with master_csv.open("r", newline="") as fh:
                        reader = csv.DictReader(fh)
                        first_row = next(reader, None)
                        if first_row:
                            panel = first_row.get("analysis_panel", "").strip()
                            if panel:
                                target_panel = panel
                except Exception:
                    pass  # Use default
            
            already_finalized = sample_id in self._finalized_samples

            # Trigger finalization (and optional SNP analysis) in background thread to avoid blocking GUI
            import threading
            def finalize_in_background():
                try:
                    logging.debug(
                        "target_bam_finalize: worker start sample_id=%s trigger_snp=%s "
                        "already_finalized=%s",
                        sample_id,
                        trigger_snp,
                        already_finalized,
                    )
                    self._set_pipeline_status(
                        sample_id,
                        phase="Finalizing target.bam",
                        progress=0.2,
                        detail="Background worker started",
                    )
                    finalization_succeeded = True

                    if already_finalized:
                        logging.info(f"Sample {sample_id} already finalized; skipping target.bam merge")
                    else:
                        logging.debug(
                            "target_bam_finalize: merge path sample_id=%s", sample_id
                        )
                        logging.info(
                            f"Triggering target.bam finalization for sample {sample_id} (status: Complete)"
                        )
                        workflow_runner = getattr(self, "workflow_runner", None)
                        if workflow_runner and hasattr(
                            workflow_runner, "submit_target_bam_finalize_job"
                        ):
                            logging.debug(
                                "target_bam_finalize: submit_target_bam_finalize_job "
                                "sample_id=%s",
                                sample_id,
                            )
                            submitted = self._resolve_submission_result(
                                workflow_runner.submit_target_bam_finalize_job(
                                    sample_dir=str(sample_dir),
                                    sample_id=sample_id,
                                    target_panel=target_panel,
                                    reference=reference_path_for_finalize,
                                ),
                                operation=f"Target BAM finalize submission for {sample_id}",
                            )
                            logging.debug(
                                "target_bam_finalize: submit finalize returned %s",
                                submitted,
                            )
                            if submitted:
                                self._set_pipeline_status(
                                    sample_id,
                                    phase="Finalize queued",
                                    progress=0.35,
                                    detail="Queued in slow worker",
                                )
                                logging.info(
                                    f"Target BAM finalization job submitted for {sample_id}"
                                )
                                self.send_update(
                                    UpdateType.WARNING_NOTIFICATION,
                                    {
                                        "title": "Target BAM finalization queued",
                                        "message": "Finalization job submitted to the workflow queue.",
                                        "sample_id": sample_id,
                                        "level": "info",
                                    },
                                    priority=6,
                                )
                                # Treat queued finalization as success; wait for target.bam before SNP submission.
                                finalization_succeeded = True
                                if trigger_snp:
                                    ready = self._wait_for_target_bam_or_timeout(
                                        sample_id, poll_s=2.0, max_wait_s=3600.0
                                    )
                                    logging.debug(
                                        "target_bam_finalize: wait target.bam sample_id=%s "
                                        "ready=%s",
                                        sample_id,
                                        ready,
                                    )
                                    if not ready:
                                        finalization_succeeded = False
                                        self._set_pipeline_status(
                                            sample_id,
                                            phase="Finalize timeout",
                                            progress=0.35,
                                            detail="target.bam not observed after queueing",
                                            level="warning",
                                        )
                                        self.send_update(
                                            UpdateType.WARNING_NOTIFICATION,
                                            {
                                                "title": "Finalization timed out",
                                                "message": "target.bam was not detected after queued finalization.",
                                                "sample_id": sample_id,
                                                "level": "warning",
                                            },
                                            priority=6,
                                        )
                            else:
                                logging.warning(
                                    f"Target BAM finalization job submission failed for {sample_id}; running inline"
                                )
                        if finalization_succeeded:
                            logging.debug(
                                "target_bam_finalize: finalize_accumulation_for_sample "
                                "sample_id=%s",
                                sample_id,
                            )
                            result = finalize_accumulation_for_sample(
                                sample_id=sample_id,
                                work_dir=str(self.monitored_directory),
                                target_panel=target_panel,
                            )
                            logging.debug(
                                "target_bam_finalize: accumulation result=%s", result
                            )
                            if result.get("status") != "error":
                                self._finalized_samples.add(sample_id)
                                self._set_pipeline_status(
                                    sample_id,
                                    phase="Finalize complete",
                                    progress=0.6 if trigger_snp else 1.0,
                                    detail="target.bam ready",
                                )
                                logging.info(
                                    f"Successfully finalized target.bam for sample {sample_id}"
                                )
                            else:
                                finalization_succeeded = False
                                self._set_pipeline_status(
                                    sample_id,
                                    phase="Finalize failed",
                                    progress=0.2,
                                    detail=str(result.get("error", "unknown error"))[:120],
                                    level="negative",
                                )
                                logging.warning(
                                    f"Target.bam finalization returned error for {sample_id}: {result.get('error', 'Unknown')}"
                                )

                    if trigger_snp and finalization_succeeded:
                        self._set_pipeline_status(
                            sample_id,
                            phase="SNP preparation",
                            progress=0.7,
                            detail="Checking prerequisites",
                        )
                        logging.debug(
                            "target_bam_finalize: SNP submission path sample_id=%s",
                            sample_id,
                        )
                        # Ensure only one SNP analysis submission runs at a time
                        with self._snp_analysis_lock:
                            running_sample = self._snp_analysis_running_sample
                            if running_sample:
                                if running_sample == sample_id:
                                    message = "SNP analysis is already running for this sample."
                                else:
                                    message = (
                                        f"SNP analysis already running for sample {running_sample}. "
                                        "Please wait for it to finish before starting another."
                                    )
                                logging.info(message)
                                self.send_update(
                                    UpdateType.WARNING_NOTIFICATION,
                                    {
                                        "title": "SNP analysis in progress",
                                        "message": message,
                                        "sample_id": sample_id,
                                        "level": "warning",
                                    },
                                    priority=6,
                                )
                                return

                            self._snp_analysis_running_sample = sample_id
                            logging.debug(
                                "target_bam_finalize: snp lock set sample_id=%s", sample_id
                            )

                        try:
                            reference_path = self._locate_reference_for_sample(sample_dir)
                            if not reference_path:
                                logging.debug(
                                    "target_bam_finalize: no reference sample_id=%s",
                                    sample_id,
                                )
                                self._set_pipeline_status(
                                    sample_id,
                                    phase="SNP skipped",
                                    progress=0.7,
                                    detail="No reference genome configured",
                                    level="warning",
                                )
                                message = (
                                    "Reference genome not available; skipping SNP analysis. "
                                    "Please provide a reference via the workflow CLI."
                                )
                                logging.warning(message)
                                self.send_update(
                                    UpdateType.WARNING_NOTIFICATION,
                                    {
                                        "title": "SNP analysis skipped",
                                        "message": message,
                                        "sample_id": sample_id,
                                        "level": "warning",
                                    },
                                    priority=6,
                                )
                                return

                            target_bam = sample_dir / "target.bam"
                            if not target_bam.exists():
                                logging.debug(
                                    "target_bam_finalize: target.bam missing sample_id=%s",
                                    sample_id,
                                )
                                self._set_pipeline_status(
                                    sample_id,
                                    phase="SNP skipped",
                                    progress=0.7,
                                    detail="target.bam missing",
                                    level="warning",
                                )
                                message = "target.bam not found; cannot run SNP analysis."
                                logging.warning(message)
                                self.send_update(
                                    UpdateType.WARNING_NOTIFICATION,
                                    {
                                        "title": "SNP analysis skipped",
                                        "message": message,
                                        "sample_id": sample_id,
                                        "level": "warning",
                                    },
                                    priority=6,
                                )
                                return

                            from robin.analysis.target_analysis import is_docker_available_for_snp_analysis
                            docker_ok, docker_error = is_docker_available_for_snp_analysis()
                            if not docker_ok:
                                logging.debug(
                                    "target_bam_finalize: docker unavailable sample_id=%s: %s",
                                    sample_id,
                                    docker_error,
                                )
                                self._set_pipeline_status(
                                    sample_id,
                                    phase="SNP skipped",
                                    progress=0.7,
                                    detail="Docker unavailable",
                                    level="warning",
                                )
                                message = (
                                    f"SNP analysis requires Docker, but it is not available. {docker_error} "
                                    "Please install Docker and ensure the daemon is running, then run SNP analysis separately."
                                )
                                logging.warning(message)
                                self.send_update(
                                    UpdateType.WARNING_NOTIFICATION,
                                    {
                                        "title": "SNP analysis skipped",
                                        "message": message,
                                        "sample_id": sample_id,
                                        "level": "warning",
                                    },
                                    priority=6,
                                )
                                return

                            workflow_runner = getattr(self, "workflow_runner", None)
                            if workflow_runner and hasattr(
                                workflow_runner, "submit_snp_analysis_job"
                            ):
                                logging.debug(
                                    "target_bam_finalize: submit_snp_analysis_job "
                                    "sample_id=%s",
                                    sample_id,
                                )
                                submitted = self._resolve_submission_result(
                                    workflow_runner.submit_snp_analysis_job(
                                        sample_dir=str(sample_dir),
                                        sample_id=sample_id,
                                        reference=reference_path,
                                        threads=4,
                                        force_regenerate=False,
                                    ),
                                    operation=f"SNP submission for {sample_id}",
                                )
                                logging.debug(
                                    "target_bam_finalize: submit_snp_analysis_job "
                                    "sample_id=%s submitted=%s",
                                    sample_id,
                                    submitted,
                                )
                            elif workflow_runner and hasattr(
                                workflow_runner, "submit_sample_job"
                            ):
                                logging.debug(
                                    "target_bam_finalize: fallback submit_sample_job "
                                    "(snp_analysis) sample_id=%s",
                                    sample_id,
                                )
                                submitted = self._resolve_submission_result(
                                    workflow_runner.submit_sample_job(
                                        sample_dir=str(sample_dir),
                                        job_type="snp_analysis",
                                        sample_id=sample_id,
                                    ),
                                    operation=f"SNP fallback submission for {sample_id}",
                                )
                                logging.debug(
                                    "target_bam_finalize: fallback submit sample_id=%s "
                                    "submitted=%s",
                                    sample_id,
                                    submitted,
                                )
                            else:
                                submitted = False

                            if submitted:
                                self._set_pipeline_status(
                                    sample_id,
                                    phase="SNP queued",
                                    progress=0.85,
                                    detail="Queued in slow worker",
                                )
                                def _mark_snp_completion(_sid: str) -> None:
                                    finished = self._wait_for_snp_outputs_or_timeout(
                                        _sid, poll_s=5.0, max_wait_s=86400.0
                                    )
                                    if finished:
                                        self._set_pipeline_status(
                                            _sid,
                                            phase="SNP complete",
                                            progress=1.0,
                                            detail="SNP outputs detected",
                                        )
                                    else:
                                        self._set_pipeline_status(
                                            _sid,
                                            phase="SNP still running",
                                            progress=0.9,
                                            detail="Waiting for outputs",
                                            level="warning",
                                        )

                                threading.Thread(
                                    target=_mark_snp_completion,
                                    args=(sample_id,),
                                    daemon=True,
                                ).start()
                                logging.info(
                                    f"SNP analysis job submitted for {sample_id} using reference {reference_path}"
                                )
                                self.send_update(
                                    UpdateType.WARNING_NOTIFICATION,
                                    {
                                        "title": "SNP analysis queued",
                                        "message": "SNP analysis job submitted to the workflow queue.",
                                        "sample_id": sample_id,
                                        "level": "info",
                                    },
                                    priority=6,
                                )
                            else:
                                self._set_pipeline_status(
                                    sample_id,
                                    phase="SNP submit failed",
                                    progress=0.75,
                                    detail="Workflow runner unavailable",
                                    level="warning",
                                )
                                logging.warning(
                                    "No workflow runner available; SNP analysis not submitted."
                                )
                                self.send_update(
                                    UpdateType.WARNING_NOTIFICATION,
                                    {
                                        "title": "SNP analysis skipped",
                                        "message": "Workflow runner unavailable; SNP analysis not submitted.",
                                        "sample_id": sample_id,
                                        "level": "warning",
                                    },
                                    priority=6,
                                )
                        except Exception as snp_exc:
                            self._set_pipeline_status(
                                sample_id,
                                phase="SNP failed",
                                progress=0.75,
                                detail=str(snp_exc)[:120],
                                level="negative",
                            )
                            logging.error(
                                f"Error submitting SNP analysis for {sample_id}: {snp_exc}",
                                exc_info=True,
                            )
                            self.send_update(
                                UpdateType.WARNING_NOTIFICATION,
                                {
                                    "title": "SNP analysis failed",
                                    "message": f"Error submitting SNP analysis: {snp_exc}",
                                    "sample_id": sample_id,
                                    "level": "negative",
                                },
                                priority=6,
                            )
                        finally:
                            with self._snp_analysis_lock:
                                self._snp_analysis_running_sample = None
                            logging.debug(
                                "target_bam_finalize: snp lock cleared sample_id=%s",
                                sample_id,
                            )
                except Exception as e:
                    logging.debug(
                        "target_bam_finalize: worker exception sample_id=%s: %s",
                        sample_id,
                        e,
                    )
                    self._set_pipeline_status(
                        sample_id,
                        phase="Finalize/SNP failed",
                        progress=0.1,
                        detail=str(e)[:120],
                        level="negative",
                    )
                    logging.error(f"Error finalizing target.bam for {sample_id}: {e}")
            
            thread = threading.Thread(target=finalize_in_background, daemon=True)
            logging.debug(
                "target_bam_finalize: starting worker thread sample_id=%s", sample_id
            )
            thread.start()

        except Exception as e:
            logging.debug(
                "target_bam_finalize: outer exception sample_id=%s: %s", sample_id, e
            )
            logging.error(f"Failed to trigger target.bam finalization for {sample_id}: {e}")
    
    def _get_expected_completion_job_types(self) -> Set[str]:
        """Get the set of job types expected to complete for this workflow."""
        if not self.workflow_steps:
            return set()
        expected: Set[str] = set()
        for step in self.workflow_steps:
            step_name = step.split(":")[-1] if ":" in step else step
            if step_name in COMPLETION_JOB_PATTERNS:
                expected.add(step_name)
        return expected

    def _expected_jobs_completed(self, sample_dir: Path, expected_job_types: Set[str]) -> bool:
        """Return True if all expected job types have completion outputs."""
        if not expected_job_types:
            return True
        if not sample_dir.exists():
            return False
        counts = self._calculate_job_counts_from_files(
            sample_dir, expected_job_types=expected_job_types
        )
        total_jobs = int(counts.get("total_jobs", 0))
        completed_jobs = int(counts.get("completed_jobs", 0))
        return total_jobs > 0 and completed_jobs >= total_jobs

    def _calculate_job_counts_from_files(
        self, sample_dir: Path, expected_job_types: Optional[Set[str]] = None
    ) -> Dict[str, Any]:
        """Calculate job counts and types from actual analysis result files in the sample directory."""
        try:
            total_jobs = 0
            completed_jobs = 0
            failed_jobs = 0
            job_types = set()
            
            job_patterns = COMPLETION_JOB_PATTERNS
            job_types_to_check = (
                expected_job_types if expected_job_types else set(job_patterns.keys())
            )
            
            # Check for each job type
            for job_type in job_types_to_check:
                patterns = job_patterns.get(job_type, [])
                job_found = False
                for pattern in patterns:
                    if "*" in pattern:
                        # Use glob pattern
                        matching_files = list(sample_dir.glob(pattern))
                        if matching_files:
                            job_found = True
                            break
                    else:
                        # Check for exact file
                        if (sample_dir / pattern).exists():
                            job_found = True
                            break
                
                if job_found:
                    total_jobs += 1
                    completed_jobs += 1  # Assume completed if files exist
                    job_types.add(job_type)

            if expected_job_types:
                total_jobs = len(expected_job_types)
            
            return {
                "total_jobs": total_jobs,
                "completed_jobs": completed_jobs,
                "failed_jobs": failed_jobs,
                "job_types": ", ".join(sorted(job_types)) if job_types else ""
            }
            
        except Exception as e:
            logging.warning(f"Error calculating job counts from files for {sample_dir}: {e}")
            return {
                "total_jobs": 0,
                "completed_jobs": 0,
                "failed_jobs": 0,
                "job_types": ""
            }

    def _create_watched_folders_page(self):
        """Create the watched folders management page."""
        try:
            from robin.workflow_ray import (
                add_watch_path,
                remove_watch_path,
                get_watched_paths,
            )
        except ImportError:
            with theme.frame(
                "R.O.B.I.N - Watched Folders",
                smalltitle="Watched Folders",
                batphone=False,
                center=self.center,
                setup_notifications=self._setup_notification_system,
            ):
                ui.notify(
                    "Manage folders requires Ray workflow. Are you running with --use-ray?",
                    type="negative",
                )
                with ui.element("div").classes("w-full min-w-0").props(
                    "id=watched-folders-page"
                ):
                    with ui.column().classes(
                        "w-full max-w-2xl mx-auto gap-3 p-2 md:p-3"
                    ):
                        with ui.element("div").classes(
                            "classification-insight-shell w-full min-w-0"
                        ):
                            ui.label("Watched folders").classes(
                                "classification-insight-heading text-headline-small"
                            )
                            with ui.element("div").classes(
                                "classification-insight-card w-full min-w-0"
                            ):
                                with ui.column().classes(
                                    "w-full min-w-0 gap-2 p-2 md:p-3"
                                ):
                                    ui.label(
                                        "Manage folders requires the Ray workflow. "
                                        "Start the app with --use-ray to enable this feature."
                                    ).classes("classification-insight-foot")
            return

        with theme.frame(
            "R.O.B.I.N - Watched Folders",
            smalltitle="Watched Folders",
            batphone=False,
            center=self.center,
            setup_notifications=self._setup_notification_system,
        ):
            with ui.element("div").classes("w-full min-w-0").props(
                "id=watched-folders-page"
            ):
                with ui.column().classes(
                    "w-full max-w-2xl mx-auto gap-3 p-2 md:p-3"
                ):
                    with ui.element("div").classes(
                        "classification-insight-shell w-full min-w-0"
                    ):
                        ui.label("Watched folders").classes(
                            "classification-insight-heading text-headline-small"
                        )
                        with ui.element("div").classes(
                            "classification-insight-card w-full min-w-0"
                        ):
                            with ui.column().classes(
                                "w-full min-w-0 gap-3 p-2 md:p-3"
                            ):
                                with ui.row().classes(
                                    "items-center gap-2 min-w-0"
                                ):
                                    ui.icon("folder_special").classes(
                                        "classification-insight-icon"
                                    )
                                    ui.label("Manage watched directories").classes(
                                        "classification-insight-model flex-1 min-w-0"
                                    )
                                ui.label(
                                    "Add or remove directories containing BAM files. "
                                    "Watched folders are monitored for new files and existing "
                                    "files are submitted for processing."
                                ).classes("classification-insight-foot")

                                watched_container = ui.column().classes(
                                    "w-full gap-2 min-w-0"
                                )
                                with watched_container:
                                    paths = get_watched_paths()
                                    if paths:
                                        ui.label("Currently watched").classes(
                                            "target-coverage-panel__meta-label mt-1 mb-1"
                                        )
                                        for p in paths:
                                            with ui.row().classes(
                                                "w-full items-center gap-2 p-2 "
                                                "rounded min-w-0 watched-folders-path-row"
                                            ):
                                                ui.label(p).classes(
                                                    "text-sm flex-1 break-all font-mono"
                                                )
                                                ui.button(
                                                    "Remove",
                                                    on_click=lambda path=p: self._do_remove_folder(
                                                        path,
                                                        watched_container,
                                                        remove_watch_path,
                                                        get_watched_paths,
                                                    ),
                                                ).props("flat color=negative size=sm no-caps")
                                    else:
                                        ui.label(
                                            "No folders currently watched."
                                        ).classes(
                                            "classification-insight-foot italic"
                                        )

                                ui.separator().classes("mgmt-detail-separator")

                                ui.label("Add folder").classes(
                                    "target-coverage-panel__meta-label mt-1 mb-1"
                                )
                                with ui.row().classes(
                                    "w-full gap-2 items-end flex-wrap"
                                ):
                                    path_input = ui.input(
                                        placeholder="e.g. /data/incoming/new_batch",
                                        label="Folder path",
                                    ).classes("flex-1 min-w-[12rem]")

                                    async def pick_folder():
                                        from robin.gui.components.folder_picker import (
                                            local_folder_picker,
                                        )

                                        watched = get_watched_paths()
                                        start = (
                                            watched[0]
                                            if watched
                                            else (
                                                self.monitored_directory
                                                or str(Path.home())
                                            )
                                        )
                                        picker = local_folder_picker(
                                            start, upper_limit=None
                                        )
                                        result = await picker
                                        if result and len(result) > 0:
                                            path_input.value = result[0]

                                    ui.button(
                                        "Browse",
                                        on_click=pick_folder,
                                        icon="folder_open",
                                    ).props("flat no-caps outline")

                                if self.monitored_directory:
                                    ui.label("Work directory").classes(
                                        "target-coverage-panel__meta-label mt-2 mb-1"
                                    )
                                    ui.label(self.monitored_directory).classes(
                                        "text-xs break-all font-mono "
                                        "watched-folders-work-dir"
                                    )

                                async def do_add_folder():
                                    path_val = (path_input.value or "").strip()
                                    if not path_val:
                                        ui.notify(
                                            "Please enter a folder path",
                                            type="warning",
                                        )
                                        return
                                    with ui.dialog().props("persistent").classes(
                                        "w-full max-w-sm"
                                    ) as add_dialog:
                                        with ui.card().classes(
                                            "robin-dialog-surface p-4 md:p-5 w-full"
                                        ):
                                            with ui.row().classes(
                                                "items-center gap-3 min-w-0"
                                            ):
                                                ui.spinner(size="lg")
                                                ui.label("Adding folder…").classes(
                                                    "classification-insight-model"
                                                )
                                    add_dialog.open()
                                    # Close modal immediately before long-running add to avoid
                                    # "client has been deleted" when user navigates away during add
                                    add_dialog.close()
                                    self._safe_notify("Adding folder...", "info")
                                    await self._do_add_folder(
                                        path_input=path_input,
                                        add_watch_path=add_watch_path,
                                        watched_container=watched_container,
                                        get_watched_paths=get_watched_paths,
                                        remove_watch_path=remove_watch_path,
                                    )

                                ui.button(
                                    "Add folder",
                                    on_click=do_add_folder,
                                    icon="add_circle_outline",
                                ).props("color=primary no-caps").classes("w-full mt-2")

    def _save_sample_identifier_manifest(
        self,
        sample_id: str,
        test_id: str,
        first_name: str,
        last_name: str,
        dob: str,
        nhs_number: str = "",
    ) -> tuple[bool, str]:
        """
        Create the sample output folder and write the identifier manifest.
        test_id is stored in plain text; first_name, last_name, dob and hospital number (nhs_number)
        are encrypted with a key derived from the date of birth.
        Returns (success, message).
        """
        if not (self.monitored_directory or "").strip():
            return False, "No output directory configured (work directory not set)."
        base = Path(self.monitored_directory)
        try:
            sample_dir = base / sample_id
            sample_dir.mkdir(parents=True, exist_ok=True)
        except OSError as e:
            return False, f"Cannot create sample folder: {e}"

        try:
            first_name_enc = _encrypt_identifier_manifest_field(first_name, dob)
            last_name_enc = _encrypt_identifier_manifest_field(last_name, dob)
            dob_enc = _encrypt_identifier_manifest_field(dob, dob)
            nhs_number_enc = _encrypt_identifier_manifest_field(nhs_number, dob) if nhs_number else ""
        except Exception as e:
            return False, f"Encryption failed: {e}"

        manifest: Dict[str, Any] = {
            "sample_id": sample_id,
            "test_id": test_id,
            "first_name": first_name_enc,
            "last_name": last_name_enc,
            "dob": dob_enc,
            "created_utc": datetime.utcnow().isoformat() + "Z",
        }
        if nhs_number_enc:
            manifest["nhs_number"] = nhs_number_enc

        manifest_path = sample_dir / SAMPLE_IDENTIFIER_MANIFEST_FILENAME
        try:
            with open(manifest_path, "w", encoding="utf-8") as f:
                json.dump(manifest, f, indent=2)
        except OSError as e:
            return False, f"Cannot write manifest: {e}"
        return True, f"Manifest saved to {sample_dir}"

    def _create_sample_id_generator_page(self):
        """Create the page for generating sample identifiers from Test ID, name, and D.O.B."""
        with theme.frame(
            "R.O.B.I.N - Generate Sample Identifier",
            smalltitle="Sample ID Generator",
            batphone=False,
            center=self.center,
            setup_notifications=self._setup_notification_system,
        ):
            with ui.element("div").classes("w-full min-w-0").props(
                "id=sample-id-generator-page"
            ):
                with ui.column().classes(
                    "w-full max-w-2xl mx-auto gap-3 p-2 md:p-3"
                ):
                    with ui.element("div").classes(
                        "classification-insight-shell w-full min-w-0"
                    ):
                        ui.label("Sample ID generator").classes(
                            "classification-insight-heading text-headline-small"
                        )
                        with ui.element("div").classes(
                            "classification-insight-card w-full min-w-0"
                        ):
                            with ui.column().classes(
                                "w-full min-w-0 gap-3 p-2 md:p-3"
                            ):
                                with ui.row().classes(
                                    "items-center gap-2 min-w-0"
                                ):
                                    ui.icon("fingerprint").classes(
                                        "classification-insight-icon"
                                    )
                                    ui.label("Generate a sample identifier").classes(
                                        "classification-insight-model flex-1 min-w-0"
                                    )
                                ui.label(
                                    "Test ID is required. First name, last name, date of birth, "
                                    "and hospital number are optional. The sample name is the MD5 "
                                    "hash of test ID, first name, last name, and D.O.B. "
                                    "(pipe-separated); empty optional fields use an empty string."
                                ).classes("classification-insight-foot")

                                test_id = ui.input(
                                    label="Test ID (required)",
                                    placeholder="e.g. LAB-2024-001",
                                ).classes("w-full")

                                first_name = ui.input(
                                    label="First name (optional)",
                                    placeholder="Given name",
                                ).classes("w-full")

                                last_name = ui.input(
                                    label="Last name (optional)",
                                    placeholder="Family name",
                                ).classes("w-full")

                                dob = ui.date_input(
                                    "Date of birth (optional)",
                                    value=None,
                                ).classes("w-full")

                                nhs_number = ui.input(
                                    label="Hospital number (optional)",
                                    placeholder="e.g. 123 456 7890",
                                ).classes("w-full")

                                ui.separator().classes("mgmt-detail-separator")

                                ui.label("Generated identifier").classes(
                                    "target-coverage-panel__meta-label mt-1 mb-1"
                                )
                                result_label = ui.label("").classes(
                                    "text-sm font-mono break-all p-3 rounded sample-id-gen-hash-preview"
                                )
                                result_label.set_visibility(False)
                                result_input = ui.input(
                                    label="Sample ID (MD5)",
                                    placeholder="Click Generate to create an ID",
                                ).classes("w-full font-mono").props("readonly outlined dense")

                                def generate_sample_id():
                                    _dob_val = dob.value
                                    dob_str = (
                                        _dob_val.strftime("%Y-%m-%d")
                                        if hasattr(_dob_val, "strftime")
                                        else (str(_dob_val).strip() if _dob_val else "")
                                    )
                                    parts = [
                                        (test_id.value or "").strip(),
                                        (first_name.value or "").strip(),
                                        (last_name.value or "").strip(),
                                        dob_str,
                                    ]
                                    if not parts[0]:
                                        ui.notify("Please enter Test ID.", type="warning")
                                        return
                                    test_id_val, first_name_val, last_name_val, dob_val = parts
                                    payload = "|".join(parts)
                                    sample_id = hashlib.md5(
                                        payload.encode("utf-8")
                                    ).hexdigest()
                                    result_input.value = sample_id
                                    result_label.set_text(f"MD5 of: {payload!r}")
                                    result_label.set_visibility(True)
                                    nhs_val = (nhs_number.value or "").strip()
                                    success, msg = self._save_sample_identifier_manifest(
                                        sample_id,
                                        test_id_val,
                                        first_name_val,
                                        last_name_val,
                                        dob_val,
                                        nhs_number=nhs_val,
                                    )
                                    if success:
                                        ui.notify(
                                            f"Sample ID generated. {msg}",
                                            type="positive",
                                        )
                                    else:
                                        ui.notify(
                                            f"Sample ID generated. {msg}",
                                            type="warning",
                                        )

                                def copy_to_clipboard():
                                    if result_input.value:
                                        ui.run_javascript(
                                            f"navigator.clipboard.writeText({json.dumps(result_input.value)})"
                                        )
                                        ui.notify("Copied to clipboard", type="positive")
                                    else:
                                        ui.notify("Generate an ID first.", type="warning")

                                with ui.row().classes(
                                    "w-full gap-2 mt-2 flex-wrap"
                                ):
                                    ui.button(
                                        "Generate sample ID",
                                        on_click=generate_sample_id,
                                        icon="fingerprint",
                                    ).props("color=primary no-caps")
                                    ui.button(
                                        "Copy to clipboard",
                                        on_click=copy_to_clipboard,
                                        icon="content_copy",
                                    ).props("flat no-caps outline")

    def _safe_notify(self, message: str, type_: str = "info"):
        """Notify only if the client context is still valid (avoids 'client deleted' errors)."""
        try:
            ui.notify(message, type=type_)
        except RuntimeError as e:
            if "deleted" not in str(e).lower():
                raise

    async def _do_add_folder(
        self,
        path_input,
        add_watch_path,
        watched_container=None,
        get_watched_paths=None,
        remove_watch_path=None,
    ):
        """Validate and add the folder path to the workflow watch (non-blocking)."""
        path_val = (path_input.value or "").strip()
        if not path_val:
            self._safe_notify("Please enter a folder path", "warning")
            return

        try:
            from nicegui import run as ng_run
            success, message = await ng_run.io_bound(add_watch_path, path_val)
        except ImportError:
            success, message = add_watch_path(path_val)

        if success:
            # If the watch add partially succeeded (e.g. some samples/subfolders were skipped),
            # surface that as a warning so users get clear visual feedback.
            try:
                msg_l = (message or "").lower()
            except Exception:
                msg_l = ""
            notify_type = (
                "warning"
                if ("skipped previously-analysed" in msg_l or "skipped previously analyzed" in msg_l)
                else "positive"
            )
            self._safe_notify(message, notify_type)
            try:
                path_input.value = ""
            except (RuntimeError, Exception) as e:
                if "deleted" not in str(e).lower() and "client" not in str(e).lower():
                    raise
                # Dialog/slot was closed or client disconnected; skip clearing input
            if watched_container is not None and get_watched_paths is not None and remove_watch_path is not None:
                try:
                    self._refresh_watched_list(watched_container, remove_watch_path, get_watched_paths)
                except RuntimeError as e:
                    if "deleted" not in str(e).lower():
                        raise
        else:
            self._safe_notify(message, "negative")

    def _do_remove_folder(self, path, watched_container, remove_watch_path, get_watched_paths):
        """Remove a folder from the watch list."""
        success, message = remove_watch_path(path)
        if success:
            ui.notify(message, type="positive")
            self._refresh_watched_list(watched_container, remove_watch_path, get_watched_paths)
        else:
            ui.notify(message, type="negative")

    def _refresh_watched_list(self, watched_container, remove_watch_path, get_watched_paths):
        """Refresh the list of watched paths in the dialog."""
        watched_container.clear()
        with watched_container:
            paths = get_watched_paths()
            if paths:
                ui.label("Currently watched").classes(
                    "target-coverage-panel__meta-label mt-1 mb-1"
                )
                for p in paths:
                    with ui.row().classes(
                        "w-full items-center gap-2 p-2 rounded min-w-0 "
                        "watched-folders-path-row"
                    ):
                        ui.label(p).classes("text-sm flex-1 break-all font-mono")
                        ui.button(
                            "Remove",
                            on_click=lambda path=p: self._do_remove_folder(
                                path,
                                watched_container,
                                remove_watch_path,
                                get_watched_paths,
                            ),
                        ).props("flat color=negative size=sm no-caps")
            else:
                ui.label("No folders currently watched.").classes(
                    "classification-insight-foot italic"
                )

    def stop_gui(self):
        """Stop the GUI thread."""
        self.is_running = False
        # Note: NiceGUI doesn't have a clean shutdown method
        # The thread will terminate when the main process ends
        logging.info("GUI shutdown requested")

    def is_gui_running(self) -> bool:
        """Check if the GUI thread is running."""
        return self.is_running and self.gui_thread and self.gui_thread.is_alive()

    def get_gui_url(self) -> str:
        """Get the URL where the GUI is running."""
        return f"http://{self.host}:{self.port}"

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
            if not hasattr(self, "workflow_runner") or self.workflow_runner is None:
                print(f"[GUI] No workflow runner available for {job_type} job")
                return False

            if sample_id is None:
                sample_id = Path(sample_dir).name

            # Check if it's a Simple workflow or Ray workflow
            if hasattr(self.workflow_runner, "submit_sample_job"):
                # Simple workflow
                return self.workflow_runner.submit_sample_job(
                    sample_dir, job_type, sample_id
                )
            elif hasattr(self.workflow_runner, "manager") and hasattr(
                self.workflow_runner.manager, "submit_sample_job"
            ):
                # Ray workflow - this is more complex, so we'll provide a fallback
                print(
                    f"[GUI] Ray workflow detected for {job_type} job - manual submission not yet supported"
                )
                return False
            else:
                print(f"[GUI] Unknown workflow type for {job_type} job")
                return False

        except Exception as e:
            print(f"[GUI] Failed to submit {job_type} job for sample {sample_id}: {e}")
            return False


def launch_gui(
    host="0.0.0.0",
    port: int = 8081,
    show: bool = False,
    workflow_runner: Any = None,
    workflow_steps: list = None,
    monitored_directory: str = "",
    center: str = None,
) -> GUILauncher:
    """Legacy launch function (kept for backward compatibility).

    Internally delegates to the refactored launcher in `robin.gui.app`.
    """
    from .gui.app import launch_gui as _launch  # type: ignore

    return _launch(
        host=host,
        port=port,
        show=show,
        reload=False,
        workflow_runner=workflow_runner,
        workflow_steps=workflow_steps,
        monitored_directory=monitored_directory,
        center=center,
    )


def get_gui_launcher() -> Optional[GUILauncher]:
    """Compatibility shim that delegates to robin.gui.app.get_gui_launcher."""
    try:
        from .gui.app import get_gui_launcher as _get  # type: ignore

        return _get()
    except Exception:
        return None


def send_gui_update(update_type: UpdateType, data: Dict[str, Any], priority: int = 0):
    """Compatibility shim that delegates to robin.gui.app.send_gui_update."""
    try:
        from .gui.app import send_gui_update as _send  # type: ignore

        _send(update_type, data, priority)
    except Exception:
        logging.info("[GUI] No GUI launcher available for update (update dropped)")


# Global reference retained for backward compatibility (now managed in gui.app)
_current_gui_launcher = None


if __name__ in {"__main__", "__mp_main__"}:
    """
    Command-line interface for direct GUI launching.

    Usage:
        # Launch GUI with default settings
        python gui_launcher.py

        # Launch GUI pointing to specific directory
        python gui_launcher.py test_out_priority_ray4

        # Launch GUI on specific port
        python gui_launcher.py test_out_priority_ray4 --port 8081

        # Launch GUI without opening browser
        python gui_launcher.py test_out_priority_ray4 --no-show

        # Launch GUI on different host
        python gui_launcher.py test_out_priority_ray4 --host 0.0.0.0

        # Launch GUI with auto-reload for development
        python gui_launcher.py test_out_priority_ray4 --reload
    """
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description="Direct GUI launcher for robin workflow monitoring",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
        Examples:
        python gui_launcher.py                                    # Launch with defaults
        python gui_launcher.py test_out_priority_ray4            # Monitor specific directory
        python gui_launcher.py test_out_priority_ray4 --port 8081 # Use different port
        python gui_launcher.py test_out_priority_ray4 --no-show   # Don't open browser
        python gui_launcher.py test_out_priority_ray4 --reload    # Enable auto-reload for development
                """,
    )

    parser.add_argument(
        "monitored_directory",
        nargs="?",
        default="",
        help="Directory containing sample output folders to monitor",
    )

    parser.add_argument(
        "--host",
        default="localhost",
        help="Host to bind the GUI to (default: localhost)",
    )

    parser.add_argument(
        "--port", type=int, default=8080, help="Port to run the GUI on (default: 8080)"
    )

    parser.add_argument(
        "--no-show", action="store_true", help="Don't automatically open the browser"
    )

    parser.add_argument(
        "--reload",
        action="store_true",
        help="Enable auto-reload for development (watches for file changes)",
    )

    args = parser.parse_args()

    # Validate monitored directory if provided
    if args.monitored_directory:
        monitored_path = Path(args.monitored_directory)
        if not monitored_path.exists():
            print(f"Error: Directory '{args.monitored_directory}' does not exist")
            sys.exit(1)
        if not monitored_path.is_dir():
            print(f"Error: '{args.monitored_directory}' is not a directory")
            sys.exit(1)
        print(f"Will monitor: {monitored_path.resolve()}")

    try:
        # Create and launch GUI directly to avoid relative import issues
        print("Creating GUI launcher...")

        # Create the GUI launcher directly
        launcher = GUILauncher(
            host=args.host,
            port=args.port,
            reload=args.reload,
        )

        # Set monitored directory if provided
        if args.monitored_directory:
            launcher.monitored_directory = str(Path(args.monitored_directory).resolve())

        print(f"Launching full GUI on http://{args.host}:{args.port}")
        print("Use Ctrl+C to stop the GUI")

        # Launch the GUI directly
        success = launcher.launch_gui(
            monitored_directory=args.monitored_directory,
            workflow_runner=None,
            workflow_steps=[],
        )

        if success:
            print("GUI launched successfully!")
            print("🌐 Open your browser to the URL above")
            print("Press Ctrl+C to stop the GUI")

            # Keep the main thread alive while GUI runs
            try:
                while launcher.is_gui_running():
                    time.sleep(1)
            except KeyboardInterrupt:
                print("\n🛑 Shutting down GUI...")
                launcher.stop_gui()
                print("GUI stopped")
        else:
            print("Failed to launch GUI")
            sys.exit(1)

    except KeyboardInterrupt:
        print("\nGUI stopped by user")
    except Exception as e:
        print(f"GUI failed to start: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)
