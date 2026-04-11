#!/usr/bin/env python3
"""
Master CSV Manager for robin

This module provides functionality to create and update master.csv files
for each sample, tracking comprehensive metadata across multiple BAM files.
"""

import os
import tempfile
import pandas as pd
from typing import Dict, Any
from dataclasses import dataclass
import logging
import time
import csv
from datetime import datetime, timezone
from pathlib import Path

# Per-sample file lock for master.csv read-modify-write
from robin.analysis.master_bed_generator import FileLock
from robin.classification_config import get_confidence_ui_tier


@dataclass
class MasterCSVManager:
    """Manages the creation and updating of master.csv files for each sample"""

    work_dir: str

    EXPORT_COLUMNS = [
        "exported_at",
        "run_summary_run_time",
        "run_summary_device",
        "run_summary_flow_cell",
        "run_summary_analysis_panel",
        "run_summary_basecall_model",
        "run_summary_sample_id",
        "run_summary_bam_passed",
        "run_summary_bam_failed",
        "run_summary_total_bases",
        "run_summary_mapped_reads",
        "run_summary_unmapped_reads",
        "run_summary_bam_batches",
        # v3: classification (see CSV_export.md): two first-seen offsets each
        "classification_sturgeon_label",
        "classification_sturgeon_first_seen_offset_sequencing",
        "classification_sturgeon_first_seen_offset_analysis",
        "classification_sturgeon_confidence_pct",
        "classification_sturgeon_probes",
        "classification_sturgeon_confidence_tier",
        "classification_nanodx_label",
        "classification_nanodx_first_seen_offset_sequencing",
        "classification_nanodx_first_seen_offset_analysis",
        "classification_nanodx_confidence_pct",
        "classification_nanodx_probes",
        "classification_nanodx_confidence_tier",
        "classification_pannanodx_label",
        "classification_pannanodx_first_seen_offset_sequencing",
        "classification_pannanodx_first_seen_offset_analysis",
        "classification_pannanodx_confidence_pct",
        "classification_pannanodx_probes",
        "classification_pannanodx_confidence_tier",
        "classification_random_forest_label",
        "classification_random_forest_first_seen_offset_sequencing",
        "classification_random_forest_first_seen_offset_analysis",
        "classification_random_forest_confidence_pct",
        "classification_random_forest_probes",
        "classification_random_forest_confidence_tier",
    ]

    _CLASSIFIER_EXPORT_SPECS = (
        ("sturgeon_scores.csv", "sturgeon", False),
        ("NanoDX_scores.csv", "nanodx", False),
        ("PanNanoDX_scores.csv", "pannanodx", False),
        ("random_forest_scores.csv", "random_forest", True),
    )

    _NA = "N/A"

    def _get_master_csv_lock_path(self, sample_id: str) -> str:
        """Path to per-sample lock file; caller must hold this lock for read-modify-write of master.csv."""
        sample_dir = os.path.join(self.work_dir, sample_id)
        lock_dir = os.path.join(sample_dir, "_locks")
        os.makedirs(lock_dir, exist_ok=True)
        return os.path.join(lock_dir, "master_csv.lock")

    def _write_csv_internal(self, data: Dict[str, Any], csv_path: str) -> None:
        """
        Write data to CSV using atomic temp-file + rename. Caller must hold the
        per-sample master_csv lock so that read-modify-write is serialized.
        """
        df = pd.DataFrame([data])
        csv_dir = os.path.dirname(csv_path)
        temp_fd, temp_path = tempfile.mkstemp(dir=csv_dir, prefix='.master_', suffix='.csv.tmp')
        try:
            with os.fdopen(temp_fd, 'w') as f:
                df.to_csv(f, index=False)
                f.flush()
                os.fsync(f.fileno())
            os.replace(temp_path, csv_path)
        except Exception:
            try:
                if os.path.exists(temp_path):
                    os.unlink(temp_path)
            except Exception:
                pass
            raise

    def update_master_csv(
        self, sample_id: str, bam_data: Dict[str, Any], bam_info: Dict[str, Any]
    ) -> None:
        """Update or create master.csv file for a sample with new BAM data. Uses per-sample lock for full read-modify-write."""
        try:
            sample_dir = os.path.join(self.work_dir, sample_id)
            os.makedirs(sample_dir, exist_ok=True)
            master_csv_path = os.path.join(sample_dir, "master.csv")
            lock_path = self._get_master_csv_lock_path(sample_id)

            with FileLock(lock_path, timeout=10.0):
                existing_data = self._load_existing_data(master_csv_path)
                if bam_info.get("state") == "pass":
                    existing_data = self._update_pass_counters(existing_data, bam_data)
                else:
                    existing_data = self._update_fail_counters(existing_data, bam_data)
                existing_data = self._update_combined_counters(existing_data, bam_data)
                existing_data = self._update_run_info(existing_data, bam_info)
                existing_data = self._update_bam_tracking(existing_data, bam_info)
                self._write_csv_internal(existing_data, master_csv_path)

            logger = logging.getLogger("robin.master_csv")
            logger.info(f"Updated master.csv for sample {sample_id}")

        except Exception as e:
            logger = logging.getLogger("robin.master_csv")
            logger.error(f"Error updating master.csv for {sample_id}: {e}")

    def _load_existing_data(self, csv_path: str) -> Dict[str, Any]:
        """
        Load existing data from CSV or create default structure.
        
        No read locking is used because:
        1. Writers use atomic write-and-rename, so readers never see partial files
        2. Worst case is reading slightly stale data (acceptable for monitoring)
        3. Keeps GUI responsive and non-blocking
        """
        # Start with default structure to ensure all fields are present
        data = self._get_default_structure()
        
        if os.path.exists(csv_path):
            try:
                # Read without locking - atomic writes ensure we never see partial data
                df = pd.read_csv(csv_path)
                if not df.empty:
                    csv_data = df.iloc[0].to_dict()
                    
                    # Merge CSV data with default structure (CSV data takes precedence)
                    for key, value in csv_data.items():
                        if value is not None and not (isinstance(value, float) and pd.isna(value)):
                            data[key] = value
                    
                    # Ensure string fields are properly converted to strings
                    # to prevent float objects from being passed to split() methods
                    string_fields = [
                        "devices", "basecall_models", "run_time", "flowcell_ids",
                        "run_info_run_time", "run_info_device", "run_info_model", 
                        "run_info_flow_cell", "samples_overview_job_types", "analysis_panel"
                    ]
                    for field in string_fields:
                        if field in data and data[field] is not None:
                            data[field] = str(data[field])
            except Exception as e:
                print(f"Warning: Error reading existing CSV: {e}")

        return data
    
    def _get_default_structure(self) -> Dict[str, Any]:
        return {
            "counter_bam_passed": 0,
            "counter_bam_failed": 0,
            "counter_mapped_count": 0,
            "counter_pass_mapped_count": 0,
            "counter_fail_mapped_count": 0,
            "counter_unmapped_count": 0,
            "counter_pass_unmapped_count": 0,
            "counter_fail_unmapped_count": 0,
            "counter_pass_bases_count": 0,
            "counter_fail_bases_count": 0,
            "counter_bases_count": 0,
            "counter_mapped_reads_num": 0,
            "counter_unmapped_reads_num": 0,
            "counter_pass_mapped_reads_num": 0,
            "counter_fail_mapped_reads_num": 0,
            "counter_pass_unmapped_reads_num": 0,
            "counter_fail_unmapped_reads_num": 0,
            "counter_mapped_bases": 0,
            "counter_unmapped_bases": 0,
            "counter_pass_mapped_bases": 0,
            "counter_fail_mapped_bases": 0,
            "counter_pass_unmapped_bases": 0,
            "counter_fail_unmapped_bases": 0,
            "devices": "",
            "basecall_models": "",
            "run_time": "",
            "flowcell_ids": "",
            "run_info_run_time": "",
            "run_info_device": "",
            "run_info_model": "",
            "run_info_flow_cell": "",
            "bam_tracking_counter": 0,
            "bam_tracking_total_files": 0,
            # Analysis panel information
            "analysis_panel": "",
            # Samples overview (persisted GUI table aggregates)
            "samples_overview_active_jobs": 0,
            "samples_overview_pending_jobs": 0,
            "samples_overview_total_jobs": 0,
            "samples_overview_completed_jobs": 0,
            "samples_overview_failed_jobs": 0,
            # Comma-separated unique job types seen for the sample
            "samples_overview_job_types": "",
            # Unix epoch seconds for last activity used by GUI
            "samples_overview_last_seen": 0.0,
        }

    def _update_pass_counters(
        self, data: Dict[str, Any], bam_data: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Update counters for pass BAM files"""
        data["counter_bam_passed"] += 1
        data["counter_pass_mapped_count"] += bam_data.get("mapped_reads", 0)
        data["counter_pass_unmapped_count"] += bam_data.get("unmapped_reads", 0)
        data["counter_pass_bases_count"] += bam_data.get("yield_tracking", 0)
        data["counter_pass_mapped_reads_num"] += bam_data.get("mapped_reads_num", 0)
        data["counter_pass_unmapped_reads_num"] += bam_data.get("unmapped_reads_num", 0)
        data["counter_pass_mapped_bases"] += bam_data.get("mapped_bases", 0)
        data["counter_pass_unmapped_bases"] += bam_data.get("unmapped_bases", 0)
        return data

    def _update_fail_counters(
        self, data: Dict[str, Any], bam_data: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Update counters for fail BAM files"""
        data["counter_bam_failed"] += 1
        data["counter_fail_mapped_count"] += bam_data.get("mapped_reads", 0)
        data["counter_fail_unmapped_count"] += bam_data.get("unmapped_reads", 0)
        data["counter_fail_bases_count"] += bam_data.get("yield_tracking", 0)
        data["counter_fail_mapped_reads_num"] += bam_data.get("mapped_reads_num", 0)
        data["counter_fail_unmapped_reads_num"] += bam_data.get("unmapped_reads_num", 0)
        data["counter_fail_mapped_bases"] += bam_data.get("mapped_bases", 0)
        data["counter_fail_unmapped_bases"] += bam_data.get("unmapped_bases", 0)
        return data

    def _update_combined_counters(
        self, data: Dict[str, Any], bam_data: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Update combined counters for all BAM files"""
        data["counter_mapped_count"] += bam_data.get("mapped_reads", 0)
        data["counter_unmapped_count"] += bam_data.get("unmapped_reads", 0)
        data["counter_bases_count"] += bam_data.get("yield_tracking", 0)
        data["counter_mapped_reads_num"] += bam_data.get("mapped_reads_num", 0)
        data["counter_unmapped_reads_num"] += bam_data.get("unmapped_reads_num", 0)
        data["counter_mapped_bases"] += bam_data.get("mapped_bases", 0)
        data["counter_unmapped_bases"] += bam_data.get("unmapped_bases", 0)
        return data

    def _update_run_info(
        self, data: Dict[str, Any], bam_info: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Update run information arrays"""
        # Update devices
        if bam_info.get("device_position"):
            device_str = str(bam_info["device_position"])
            devices_str = str(data.get("devices", ""))
            devices = devices_str.split(",") if devices_str else []
            if device_str not in devices:
                devices.append(device_str)
            data["devices"] = ",".join(devices)
            data["run_info_device"] = device_str

        # Update basecall models
        if bam_info.get("basecall_model"):
            model_str = str(bam_info["basecall_model"])
            models_str = str(data.get("basecall_models", ""))
            models = models_str.split(",") if models_str else []
            if model_str not in models:
                models.append(model_str)
            data["basecall_models"] = ",".join(models)
            data["run_info_model"] = model_str

        # Update flow cell IDs
        if bam_info.get("flow_cell_id"):
            flowcell_str = str(bam_info["flow_cell_id"])
            flowcells_str = str(data.get("flowcell_ids", ""))
            flowcells = flowcells_str.split(",") if flowcells_str else []
            if flowcell_str not in flowcells:
                flowcells.append(flowcell_str)
            data["flowcell_ids"] = ",".join(flowcells)
            data["run_info_flow_cell"] = flowcell_str

        # Update run time
        if bam_info.get("time_of_run"):
            # Convert to string if it's a float or other type
            time_str = str(bam_info["time_of_run"])
            run_times_str = str(data.get("run_time", ""))
            run_times = run_times_str.split(",") if run_times_str else []
            if time_str not in run_times:
                run_times.append(time_str)
            data["run_time"] = ",".join(run_times)
            data["run_info_run_time"] = time_str

        return data

    def _update_bam_tracking(
        self, data: Dict[str, Any], bam_info: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Update BAM tracking information"""
        data["bam_tracking_counter"] += 1
        data["bam_tracking_total_files"] += 1

        return data

    def update_analysis_panel(self, sample_id: str, panel: str) -> None:
        """
        Update the analysis panel information for a sample.
        Uses per-sample lock for full read-modify-write.
        """
        try:
            sample_dir = os.path.join(self.work_dir, sample_id)
            os.makedirs(sample_dir, exist_ok=True)
            master_csv_path = os.path.join(sample_dir, "master.csv")
            lock_path = self._get_master_csv_lock_path(sample_id)

            with FileLock(lock_path, timeout=10.0):
                existing_data = self._load_existing_data(master_csv_path)
                existing_data["analysis_panel"] = str(panel)
                self._write_csv_internal(existing_data, master_csv_path)

            logger = logging.getLogger("robin.master_csv")
            logger.info(f"Updated analysis panel to '{panel}' for sample {sample_id}")

        except Exception as e:
            logger = logging.getLogger("robin.master_csv")
            logger.error(f"Error updating analysis panel for {sample_id}: {e}")

    # ---------------------------------------------------------------------
    # Public helpers for GUI/workflow to persist overview stats
    # ---------------------------------------------------------------------
    def update_sample_overview(self, sample_id: str, overview: Dict[str, Any]) -> None:
        """
        Persist high-level per-sample overview stats used by the GUI table.

        Expected keys in ``overview`` (all optional, defaults applied when missing):
          - active_jobs: int
          - pending_jobs: int
          - total_jobs: int
          - completed_jobs: int
          - failed_jobs: int
          - job_types: Iterable[str] | str (stored as comma-separated unique list)
          - last_seen: float | int (unix epoch seconds)
        """
        try:
            sample_dir = os.path.join(self.work_dir, sample_id)
            os.makedirs(sample_dir, exist_ok=True)
            master_csv_path = os.path.join(sample_dir, "master.csv")
            lock_path = self._get_master_csv_lock_path(sample_id)

            with FileLock(lock_path, timeout=10.0):
                existing_data = self._load_existing_data(master_csv_path)
                existing_data["samples_overview_active_jobs"] = int(
                    overview.get(
                        "active_jobs", existing_data.get("samples_overview_active_jobs", 0)
                    )
                )
                existing_data["samples_overview_pending_jobs"] = int(
                    overview.get(
                        "pending_jobs",
                        existing_data.get("samples_overview_pending_jobs", 0),
                    )
                )
                existing_data["samples_overview_total_jobs"] = int(
                    overview.get(
                        "total_jobs", existing_data.get("samples_overview_total_jobs", 0)
                    )
                )
                existing_data["samples_overview_completed_jobs"] = int(
                    overview.get(
                        "completed_jobs",
                        existing_data.get("samples_overview_completed_jobs", 0),
                    )
                )
                existing_data["samples_overview_failed_jobs"] = int(
                    overview.get(
                        "failed_jobs", existing_data.get("samples_overview_failed_jobs", 0)
                    )
                )
                jt_value = overview.get("job_types")
                if isinstance(jt_value, str):
                    new_types = {t.strip() for t in jt_value.split(",") if t.strip()}
                elif jt_value is None:
                    new_types = set()
                else:
                    try:
                        new_types = {str(t).strip() for t in jt_value if str(t).strip()}
                    except Exception:
                        new_types = set()
                if existing_data.get("samples_overview_job_types"):
                    old_types = {
                        t.strip()
                        for t in str(
                            existing_data.get("samples_overview_job_types", "")
                        ).split(",")
                        if t.strip()
                    }
                else:
                    old_types = set()
                all_types = sorted(old_types.union(new_types))
                existing_data["samples_overview_job_types"] = ",".join(all_types)
                try:
                    last_seen = float(overview.get("last_seen"))
                except Exception:
                    last_seen = float(existing_data.get("samples_overview_last_seen", 0.0))
                existing_data["samples_overview_last_seen"] = last_seen
                self._write_csv_internal(existing_data, master_csv_path)

            logger = logging.getLogger("robin.master_csv")
            logger.debug(
                f"Updated samples overview in master.csv for sample {sample_id}"
            )
        except Exception as e:
            logger = logging.getLogger("robin.master_csv")
            logger.warning(f"Failed to update samples overview for {sample_id}: {e}")

    # ---------------------------------------------------------------------
    # Master export CSV helpers
    # ---------------------------------------------------------------------
    @staticmethod
    def _get_case_insensitive(row: Dict[str, Any], key: str) -> str:
        """Get a CSV row value by key, case-insensitive."""
        for k, v in row.items():
            if str(k).strip().lower() == key.lower():
                return "" if v is None else str(v).strip()
        return ""

    @staticmethod
    def _int_string(value: Any) -> str:
        """Return integer as plain string, without separators/scientific notation."""
        if value is None:
            return ""
        s = str(value).strip().replace(",", "")
        if not s:
            return ""
        try:
            if "e" in s.lower() or "." in s:
                return str(int(float(s)))
            return str(int(s))
        except (ValueError, TypeError, OverflowError):
            return ""

    @staticmethod
    def _format_exported_at_utc() -> str:
        """UTC ISO-8601 timestamp for exported_at column."""
        return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    @staticmethod
    def _format_offset_hhmmss(delta_seconds: float) -> str:
        """Non-negative duration as H:MM:SS (hours unbounded)."""
        total = int(max(0.0, delta_seconds))
        h = total // 3600
        m = (total % 3600) // 60
        s = total % 60
        return f"{h}:{m:02d}:{s:02d}"

    @staticmethod
    def _parse_run_start_epoch(master_row: Dict[str, Any]) -> float | None:
        """Parse run_info_run_time (sequencing / BAM acquisition) to Unix epoch seconds."""
        raw = MasterCSVManager._get_case_insensitive(
            master_row, "run_info_run_time"
        )
        if not raw or not str(raw).strip():
            return None
        try:
            raw_ts = str(raw).strip().replace("Z", "+00:00")
            dt = datetime.fromisoformat(raw_ts)
            if dt.tzinfo is not None:
                dt = dt.astimezone()
            return dt.timestamp()
        except Exception:
            return None

    @staticmethod
    def _timestamp_ms_from_row(row: Dict[str, Any]) -> float | None:
        for k, v in row.items():
            if str(k).strip().lower() == "timestamp":
                try:
                    return float(v)
                except (TypeError, ValueError):
                    return None
        return None

    @staticmethod
    def _number_probes_from_row(row: Dict[str, Any]) -> int:
        for k, v in row.items():
            if str(k).strip().lower() == "number_probes":
                try:
                    if v is None or str(v).strip() == "":
                        return 0
                    return int(float(v))
                except (TypeError, ValueError):
                    return 0
        return 0

    @staticmethod
    def _is_score_column(col_name: str) -> bool:
        c = str(col_name).strip().lower()
        return c not in ("timestamp", "number_probes")

    @classmethod
    def _winner_from_row(cls, row: Dict[str, Any]) -> tuple[str, float]:
        """Argmax over class score columns; same rule as summary._extract_classification_data."""
        max_score = 0.0
        best_class = "Unknown"
        for col, value in row.items():
            if not cls._is_score_column(col):
                continue
            try:
                score = float(value)
                if score > max_score:
                    max_score = score
                    best_class = str(col).strip()
            except (TypeError, ValueError):
                continue
        return best_class, max_score

    @staticmethod
    def _normalize_confidence_pct(raw_max: float, is_random_forest: bool) -> float:
        if is_random_forest:
            return raw_max * 100 if raw_max <= 1 else raw_max
        return raw_max * 100

    @staticmethod
    def _tier_export_string(classifier_key: str, confidence_percent: float) -> str:
        t = get_confidence_ui_tier(classifier_key, float(confidence_percent))
        return {"high": "High", "medium": "Medium", "low": "Low"}.get(
            str(t).lower(), "Low"
        )

    @staticmethod
    def _format_confidence_pct_value(confidence_percent: float) -> str:
        v = float(confidence_percent)
        if abs(v - round(v)) < 1e-9:
            return str(int(round(v)))
        return f"{v:.1f}"

    @classmethod
    def _sorted_score_rows(cls, rows: list[Dict[str, Any]]) -> list[Dict[str, Any]]:
        """Rows with timestamps first (ascending by time, then file order); rows without ts after."""
        enriched: list[tuple[bool, float, int, Dict[str, Any]]] = []
        for i, r in enumerate(rows):
            ts = cls._timestamp_ms_from_row(r)
            if ts is not None:
                enriched.append((False, ts, i, r))
            else:
                enriched.append((True, 0.0, i, r))
        enriched.sort(key=lambda x: (x[0], x[1], x[2]))
        return [t[3] for t in enriched]

    @classmethod
    def _analysis_start_epoch(cls, rows: list[Dict[str, Any]]) -> float | None:
        """Earliest parseable score-row timestamp (ms → s); analysis wall-clock start."""
        best: float | None = None
        for r in rows:
            ts = cls._timestamp_ms_from_row(r)
            if ts is None:
                continue
            sec = float(ts) / 1000.0
            if best is None or sec < best:
                best = sec
        return best

    @classmethod
    def _read_score_csv_rows(cls, path: Path) -> list[Dict[str, Any]]:
        if not path.exists():
            return []
        try:
            with path.open("r", newline="", encoding="utf-8", errors="replace") as f:
                reader = csv.DictReader(f)
                return [row for row in reader if row]
        except Exception:
            return []

    def _set_classifier_na(self, out: Dict[str, str], prefix: str) -> None:
        out[f"{prefix}_label"] = self._NA
        out[f"{prefix}_first_seen_offset_sequencing"] = self._NA
        out[f"{prefix}_first_seen_offset_analysis"] = self._NA
        out[f"{prefix}_confidence_pct"] = self._NA
        out[f"{prefix}_probes"] = self._NA
        out[f"{prefix}_confidence_tier"] = self._NA

    def _merge_classification_export(
        self,
        out: Dict[str, str],
        sample_id: str,
        master_row: Dict[str, Any],
    ) -> None:
        """Populate v3 classification columns per CSV_export.md."""
        sample_dir = Path(self.work_dir) / sample_id
        t_sequencing_epoch = self._parse_run_start_epoch(master_row)

        for fname, cfg_key, is_rf in self._CLASSIFIER_EXPORT_SPECS:
            prefix = f"classification_{cfg_key}"
            score_path = sample_dir / fname
            rows = self._read_score_csv_rows(score_path)
            if not rows:
                self._set_classifier_na(out, prefix)
                continue

            last = rows[-1]
            label, raw_max = self._winner_from_row(last)
            probes = self._number_probes_from_row(last)
            conf_pct = self._normalize_confidence_pct(raw_max, is_rf)
            tier = self._tier_export_string(cfg_key, conf_pct)

            first_epoch: float | None = None
            for r in self._sorted_score_rows(rows):
                wl, _ = self._winner_from_row(r)
                if wl != label:
                    continue
                ts_ms = self._timestamp_ms_from_row(r)
                if ts_ms is None:
                    continue
                first_epoch = float(ts_ms) / 1000.0
                break

            offset_seq = self._NA
            offset_ana = self._NA
            if first_epoch is not None:
                if t_sequencing_epoch is not None:
                    offset_seq = self._format_offset_hhmmss(
                        max(0.0, first_epoch - t_sequencing_epoch)
                    )
                t_analysis = self._analysis_start_epoch(rows)
                if t_analysis is not None:
                    offset_ana = self._format_offset_hhmmss(
                        max(0.0, first_epoch - t_analysis)
                    )

            out[f"{prefix}_label"] = label
            out[f"{prefix}_first_seen_offset_sequencing"] = offset_seq
            out[f"{prefix}_first_seen_offset_analysis"] = offset_ana
            out[f"{prefix}_confidence_pct"] = self._format_confidence_pct_value(conf_pct)
            out[f"{prefix}_probes"] = str(int(probes))
            out[f"{prefix}_confidence_tier"] = tier

    def _load_master_csv_first_row(self, sample_id: str) -> Dict[str, Any]:
        """Load first row from per-sample master.csv."""
        sample_dir = os.path.join(self.work_dir, sample_id)
        master_csv_path = os.path.join(sample_dir, "master.csv")
        if not os.path.exists(master_csv_path):
            return {}
        try:
            with open(master_csv_path, "r", newline="") as f:
                reader = csv.DictReader(f)
                return next(reader, {}) or {}
        except Exception:
            return {}

    def build_master_export_row(self, sample_id: str, exported_at: str) -> Dict[str, str]:
        """Build one export row for a sample using CSV_export.md (run summary + v2 classification)."""
        row = self._load_master_csv_first_row(sample_id)
        out: Dict[str, str] = {c: "" for c in self.EXPORT_COLUMNS}

        out["exported_at"] = exported_at
        out["run_summary_sample_id"] = str(sample_id).strip()
        out["run_summary_run_time"] = self._get_case_insensitive(row, "run_info_run_time")
        out["run_summary_device"] = (
            self._get_case_insensitive(row, "run_info_device")
            or self._get_case_insensitive(row, "devices")
        )
        out["run_summary_flow_cell"] = (
            self._get_case_insensitive(row, "run_info_flow_cell")
            or self._get_case_insensitive(row, "flowcell_ids")
        )
        out["run_summary_analysis_panel"] = self._get_case_insensitive(row, "analysis_panel")
        out["run_summary_basecall_model"] = (
            self._get_case_insensitive(row, "run_info_model")
            or self._get_case_insensitive(row, "basecall_models")
        )
        out["run_summary_bam_passed"] = self._int_string(
            self._get_case_insensitive(row, "counter_bam_passed")
        )
        out["run_summary_bam_failed"] = self._int_string(
            self._get_case_insensitive(row, "counter_bam_failed")
        )
        out["run_summary_total_bases"] = self._int_string(
            self._get_case_insensitive(row, "counter_bases_count")
        )
        out["run_summary_mapped_reads"] = self._int_string(
            self._get_case_insensitive(row, "counter_mapped_count")
        )
        out["run_summary_unmapped_reads"] = self._int_string(
            self._get_case_insensitive(row, "counter_unmapped_count")
        )

        btc = self._int_string(self._get_case_insensitive(row, "bam_tracking_counter"))
        btt = self._int_string(self._get_case_insensitive(row, "bam_tracking_total_files"))
        if btc and btt:
            out["run_summary_bam_batches"] = f"{btc} / {btt}"
        elif btc:
            out["run_summary_bam_batches"] = btc
        elif btt:
            out["run_summary_bam_batches"] = btt

        self._merge_classification_export(out, sample_id, row)

        return out

    def export_master_csv_rows(self, sample_ids: list[str], output_path: str) -> str:
        """
        Export one-row-per-sample master export CSV.
        Returns output_path for convenience.
        """
        out_dir = os.path.dirname(output_path)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        exported_at = self._format_exported_at_utc()
        rows = [self.build_master_export_row(sid, exported_at) for sid in sample_ids]
        with open(output_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=self.EXPORT_COLUMNS, extrasaction="ignore")
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
        return output_path
